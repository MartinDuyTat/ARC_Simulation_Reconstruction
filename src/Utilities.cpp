// Martin Duy Tat 6th September 2022

#include"TRandom.h"
#include"TMath.h"
#include"Utilities.h"
#include"Settings.h"
#include"RadiatorCell.h"
#include"BarrelRadiatorCell.h"
#include"ParticleTrack.h"
#include"RadiatorArray.h"
#include"PhotonMapper.h"
#include"PhotonReconstructor.h"

namespace Utilities {

  Vector VectorFromSpherical(double R, double CosTheta, double Phi) {
    const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
    const double CosPhi = TMath::Cos(Phi);
    const double SinPhi = TMath::Sin(Phi);
    return Vector{R*CosPhi*SinTheta, R*SinPhi*SinTheta, R*CosTheta};
  }

  Vector GenerateRandomBarrelTrackZRange(double z_min, double z_max) {
    const double Radius = Settings::GetDouble("ARCGeometry/Radius");
    const double z = gRandom->Uniform(z_min, z_max);
    const double CosTheta = z/TMath::Sqrt(z*z + Radius*Radius);
    const double Phi = Settings::GetBool("Particle/RandomPhi")
      ? gRandom->Uniform(-TMath::Pi(), TMath::Pi())
      : gRandom->Uniform(Settings::GetDouble("Particle/Phi_min"),
			 Settings::GetDouble("Particle/Phi_max"));
    const double MomentumMag = Settings::GetDouble("Particle/Momentum");
    const Vector Momentum = VectorFromSpherical(MomentumMag, CosTheta, Phi);
    return Momentum;
  }

  Vector GenerateRandomBarrelTrack(double &CosTheta, double &Phi) {
    const double Radius = Settings::GetDouble("ARCGeometry/Radius");
    const double z = gRandom->Uniform(Settings::GetDouble("Particle/z_min"),
				      Settings::GetDouble("Particle/z_max"));
    CosTheta = z/TMath::Sqrt(z*z + Radius*Radius);
    Phi = Settings::GetBool("Particle/RandomPhi")
      ? gRandom->Uniform(-TMath::Pi(), TMath::Pi())
      : gRandom->Uniform(Settings::GetDouble("Particle/Phi_min"),
			 Settings::GetDouble("Particle/Phi_max"));
    const double MomentumMag = Settings::GetDouble("Particle/Momentum");
    const Vector Momentum = VectorFromSpherical(MomentumMag, CosTheta, Phi);
    return Momentum;
  }

  Vector GenerateRandomEndCapTrack() {
    const double Radius = Settings::GetDouble("ARCGeometry/Radius") + 0.10;
    const double z = Settings::GetDouble("ARCGeometry/BarrelZ");
    auto GetUniformCircle = [Radius] () {
      double x = Radius, y = Radius;
      while(x*x + y*y > Radius*Radius) {
	x = gRandom->Uniform(-Radius, Radius);
	y = gRandom->Uniform(-Radius, Radius);
      }
      return std::make_pair(x, y);
    };
    auto [x, y] = GetUniformCircle();
    const double MomentumMag = Settings::GetDouble("Particle/Momentum");
    const Vector Momentum = Vector(x, y, z).Unit()*MomentumMag;
    return Momentum;
  }

  ResolutionStruct TrackPhotons(ParticleTrack particleTrack,
				const RadiatorCell &radiatorCell,
				const RadiatorArray &radiatorArray) {
    if(!particleTrack.FindRadiator(radiatorArray)) {
      return ResolutionStruct{};
    }
    particleTrack.TrackThroughAerogel();
    particleTrack.TrackThroughGasToMirror();
    if(particleTrack.GetParticleLocation() == ParticleTrack::Location::MissedMirror) {
      return ResolutionStruct{0.0, 0, true, true};
    }
    std::size_t Counter = 0;
    while(particleTrack.GetParticleLocation() != ParticleTrack::Location::Mirror) {
      bool IsEdgeCell = particleTrack.GetRadiatorCell()->IsEdgeCell();
      particleTrack.ConvertBackToGlobalCoordinates();
      if(!particleTrack.FindRadiator(radiatorArray) || Counter > 10) {
	if(IsEdgeCell && radiatorCell.IsEdgeCell()) {
	  return ResolutionStruct{0.0, 0, true, true};
	} else {
	  return ResolutionStruct{};
	}
      }
      particleTrack.TrackThroughGasToMirror();
      if(particleTrack.GetParticleLocation() == ParticleTrack::Location::MissedMirror) {
	return ResolutionStruct{0.0, 0, true, true};
      }
      Counter++;
    }
    if(particleTrack.GetRadiatorCell() != &radiatorCell) {
      return ResolutionStruct{};
    }
    auto Photons = particleTrack.GetParticleLocation() != ParticleTrack::Location::Mirror ?
    std::vector<Photon>() : particleTrack.GeneratePhotonsFromGas();
    ResolutionStruct resolutionStruct;
    resolutionStruct.HitCorrectCell = true;
    std::vector<double> CherenkovAngles;
    CherenkovAngles.reserve(Photons.size());
    for(auto &Photon : Photons) {
      auto photonHit = PhotonMapper::TracePhoton(Photon, radiatorArray);
      if(Photon.GetStatus() == Photon::Status::DetectorHit) {
	auto reconstructedPhoton =
	  PhotonReconstructor::ReconstructPhoton(*photonHit);
	const double CherenkovAngle = TMath::ACos(reconstructedPhoton.m_CosCherenkovAngle);
	if(CherenkovAngle < 0.0) {
	  continue;
	}
	CherenkovAngles.push_back(CherenkovAngle);
	resolutionStruct.N++;
	resolutionStruct.CentreHitDistance += photonHit->m_CentreHitDistance;
      } else if(Photon.GetStatus() == Photon::Status::WallMiss ||
		Photon.GetStatus() == Photon::Status::Backwards ||
		Photon.GetStatus() == Photon::Status::DetectorMiss) {
	resolutionStruct.HitTopWall = true;
      }
    }
    if(resolutionStruct.N > 0) {
      resolutionStruct.x = TMath::RMS(CherenkovAngles.begin(),
				      CherenkovAngles.end());
      resolutionStruct.x /= TMath::Sqrt(resolutionStruct.N);
      resolutionStruct.CentreHitDistance /= resolutionStruct.N;
    }
    return resolutionStruct;
  }

}
