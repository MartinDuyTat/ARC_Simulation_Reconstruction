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
#include"ParticleMass.h"

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
    const Vector Momentum = VectorFromSpherical(GetMomentumMag(), CosTheta, Phi);
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
    const Vector Momentum = VectorFromSpherical(GetMomentumMag(), CosTheta, Phi);
    return Momentum;
  }

  Vector GenerateRandomEndCapTrack() {
    const double Radius = Settings::GetDouble("ARCGeometry/MaxEndCapRadius");
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
    const Vector Momentum = Vector(x, y, z).Unit()*GetMomentumMag();
    return Momentum;
  }

double GetMomentumMag() {
  if(Settings::GetBool("Particle/ConstantMomentum")) {
    return Settings::GetDouble("Particle/Momentum");
  } else {
    const double Momentum_min = Settings::GetDouble("Particle/Momentum_min");
    const double Momentum_max = Settings::GetDouble("Particle/Momentum_max");
    const double Momentum_log_min = TMath::Log(Momentum_min);
    const double Momentum_log_max = TMath::Log(Momentum_max);
    const double Momentum_log = gRandom->Uniform(Momentum_log_min,
						 Momentum_log_max);
    return TMath::Exp(Momentum_log);
  }
}

  ResolutionStruct TrackPhotons(ParticleTrack particleTrack,
				const RadiatorCell &radiatorCell,
				const RadiatorArray &radiatorArray) {
    // First find the correct cell
    if(!particleTrack.FindRadiator(radiatorArray)) {
      return ResolutionStruct{};
    }
    // Track through the cell
    particleTrack.TrackThroughAerogel();
    particleTrack.TrackThroughGasToMirror();
    if(particleTrack.GetParticleLocation() != ParticleTrack::Location::Mirror) {
      particleTrack.TrackToNextCell(radiatorArray);
    }
    // Check if track hit the cell we're optimising
    if(particleTrack.GetRadiatorCell() != &radiatorCell) {
      return ResolutionStruct{};
    }
    // If we're in the cell we're optimising, make sure we hit the mirror
    const auto Location = particleTrack.GetParticleLocation();
    if(Location == ParticleTrack::Location::MissedMirror) {
      return ResolutionStruct{0.0, 0, true, true};
    }
    // If particle hit the mirror, generate photons from gas
    auto Photons = Location == ParticleTrack::Location::Mirror ?
                   particleTrack.GeneratePhotonsFromGas() :
                   std::vector<Photon>();
    // Loop over all photons
    ResolutionStruct resolutionStruct;
    resolutionStruct.HitCorrectCell = true;
    std::vector<double> CherenkovAngles;
    CherenkovAngles.reserve(Photons.size());
    for(auto &Photon : Photons) {
      // Trace photons through optics and obtain detector hits
      auto photonHit = PhotonMapper::TracePhoton(Photon, radiatorArray);
      if(Photon.GetRadiatorCell() != &radiatorCell) {
	// If photon hits a different cell, skip
	continue;
      } else if(Photon.GetStatus() == Photon::Status::DetectorHit) {
	// If photon is detected, reconstruct Cherenkov angle
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
	// If photon doesn't hit detector, a penalty will be added later
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

  double GetIndexRefraction(Photon::Radiator Radiator,
			    bool ChromaticDispersion,
			    double Energy) {
    switch(Radiator) {
    case Photon::Radiator::Aerogel:
      {
	auto GetIndex = [] (double eph) {
	  const double InsideSqrt = 46.41/(113.8 - eph*eph) +
	  228.7/(328.5 - eph*eph);
	  return 1.0 + 0.03*2.1467*(TMath::Sqrt(1.0 + InsideSqrt) - 1.0);
	};
	if(ChromaticDispersion) {
	  // Equation from Roger Forty via email
	  return GetIndex(Energy);
	} else {
	  const double eph = 1239.841987427/400.0;
	  return GetIndex(eph);
	}
      }
    case Photon::Radiator::Gas:
      {
	auto GetIndex = [] (double L) {
	  return 1.0 + 0.25324*1e-6/((1.0/(73.7*73.7)) - (1.0/(L*L)));
	};
	if(ChromaticDispersion) {
	  // Pressure at 1.0 bar
	  // Sellmeier equation with coefficients from
	  // https://twiki.cern.ch/twiki/bin/view/LHCb/C4F10
	  // They are similar to A. Filippas, et al. Nucl. Instr. and Meth. B, 196 (2002),
	  // p. 340 but now quite...?
	  const double Lambda = 1239.841987427/Energy;
	  return GetIndex(Lambda);
	} else {
	  return GetIndex(400);
	}
      }
    default:
      return 1.0;
    }
  }

  double GetPredictedCherenkovAngle(double Momentum,
				    int ID,
				    Photon::Radiator Radiator) {
    const double n = GetIndexRefraction(Radiator, false);
    const double Mass = ParticleMass::GetMass(ID);
    const double Energy = TMath::Sqrt(Mass*Mass + Momentum*Momentum);
    const double v = Momentum/Energy;
    return TMath::ACos(1.0/(v*n));
  }

  double GetCherenkovAngleDifference(double Momentum,
				     int ID1,
				     int ID2,
				     Photon::Radiator Radiator) {
    const double Angle1 = GetPredictedCherenkovAngle(Momentum, ID1, Radiator);
    const double Angle2 = GetPredictedCherenkovAngle(Momentum, ID2, Radiator);
    return TMath::Abs(Angle1 - Angle2);
  }

}
