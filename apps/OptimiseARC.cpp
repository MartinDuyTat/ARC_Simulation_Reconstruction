// Martin Duy Tat 1st May 2022
/**
 * OptimiseARC is an application for optimising an ARC cell
 * Run with RunARC ColumnNumber RowNumber <settings name> <settings filename> ...
 * There can be an arbitrary number of settings files added
 */

#include<iostream>
#include<string>
#include"TRandom.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"
#include"Settings.h"
#include"RadiatorCell.h"
#include"TrackingVolume.h"
#include"ParticleTrack.h"
#include"PhotonMapper.h"
#include"PhotonReconstructor.h"

using Vector = ROOT::Math::XYZVector;

Vector VectorFromSpherical(double R, double CosTheta, double Phi);

double CalculateResolution(const TrackingVolume &InnerTracker,
			   RadiatorCell &radiatorCell,
			   double z_min,
			   double z_max,
			   double Radius,
			   double MomentumMag,
			   double ParticleID);
			   

int main(int argc, char *argv[]) {
  if(argc%2 != 1 || argc < 3) {
    return 0;
  }
  std::cout << "Welcome to the ARC optimiser\n";
  for(int i = 3; i < argc; i += 2) {
    const std::string SettingsName = argv[i];
    const std::string SettingsFilename = argv[i + 1];
    Settings::AddSettings(SettingsName, SettingsFilename);
    std::cout << "Added " << SettingsName << " settings from " << SettingsFilename << "\n";
  }
  const std::string RunMode(argv[1]);
  if(Settings::Exists("General/Seed")) {
    gRandom->SetSeed(Settings::GetInt("General/Seed"));
  }
  const int ParticleID = Settings::GetInt("Particle/ID");;
  const TrackingVolume InnerTracker;
  const int CellsPerRow = 2*Settings::GetDouble("ARCGeometry/CellsPerRow") - 1;
  const double HexagonSize = Settings::GetDouble("ARCGeometry/Length")/2*CellsPerRow;
  RadiatorCell radiatorCell(std::stoi(std::string(argv[1])),
			    std::stoi(std::string(argv[2])), 
			    HexagonSize);
  const Vector CellPosition = radiatorCell.GetRadiatorPosition();
  const double Radius = Settings::GetDouble("ARCGeometry/Radius");
  const double MomentumMag = Settings::GetDouble("Particle/Momentum");
  const double z_min = CellPosition.X() - HexagonSize/2.0;
  const double z_max = CellPosition.X() + HexagonSize/2.0;

  auto MinimiseFunction = [&] (const double *x) {
    radiatorCell.SetMirrorCurvature(x[0]);
    radiatorCell.SetMirrorPosition(x[1], x[2]);
    const double Resolution = CalculateResolution(InnerTracker,
						  radiatorCell,
						  z_min,
						  z_max,
						  Radius,
						  MomentumMag,
						  ParticleID);
    std::cout << Resolution << "\n";
    return Resolution;
  };
  ROOT::Math::Functor fcn(MinimiseFunction, 3);
  ROOT::Minuit2::Minuit2Minimizer Minimiser;
  Minimiser.SetPrintLevel(4);
  Minimiser.SetPrecision(0.01);
  Minimiser.SetFunction(fcn);
  Minimiser.SetVariable(0, "MirrorCurvature", 0.37, 0.01);
  Minimiser.SetVariableLimits(0, 0.25, 0.45);
  Minimiser.SetVariable(1, "Mirror_xPosition", 0.0, 0.001);
  Minimiser.SetVariableLimits(1, -0.01, 0.01);
  Minimiser.SetVariable(1, "Mirror_zPosition", 0.0, 0.001);
  Minimiser.SetVariableLimits(2, -0.01, 0.01);
  std::cout << "Starting minimisation...\n";
  Minimiser.Minimize();
  std::cout << "ARC is optimised!\n";
  return 0;
}

Vector VectorFromSpherical(double R, double CosTheta, double Phi) {
  const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
  const double CosPhi = TMath::Cos(Phi);
  const double SinPhi = TMath::Sin(Phi);
  return Vector{R*CosPhi*SinTheta, R*SinPhi*SinTheta, R*CosTheta};
}

double CalculateResolution(const TrackingVolume &InnerTracker,
			   RadiatorCell &radiatorCell,
			   double z_min,
			   double z_max,
			   double Radius,
			   double MomentumMag, 
			   double ParticleID) {
  int NumberPhotons = 0;
  // x is the Cherenkov angle, x2 is the Cherenkov angle squared
  double x = 0.0, x2 = 0.0;
  double Resolution = 0.0, NewResolution = 0.0;
  constexpr double Precision = 0.01;
  while(true) {
    for(int i = 0; i < 1000; i++) {
      const double z = gRandom->Uniform(z_min, z_max);
      const double CosTheta = z/TMath::Sqrt(z*z + Radius*Radius);
      const double Phi = gRandom->Uniform(-0.1, 0.1);
      const Vector Momentum = VectorFromSpherical(MomentumMag, CosTheta, Phi);
      ParticleTrack particleTrack(ParticleID, Momentum, Vector(0.0, 0.0, 0.0));
      particleTrack.TrackThroughTracker(InnerTracker);
      particleTrack.SetRadiator(radiatorCell);
      particleTrack.ConvertToRadiatorCoordinates();
      if(particleTrack.GetParticleLocation() != ParticleTrack::Location::EntranceWindow) {
	continue;
      }
      particleTrack.TrackThroughRadiatorCell();
      if(particleTrack.GetParticleLocation() != ParticleTrack::Location::Mirror) {
	continue;
      }
      auto Photons = particleTrack.GeneratePhotonsFromGas();
      for(auto &Photon : Photons) {
	PhotonMapper::TracePhoton(Photon);
	if(Photon.m_Status == Photon::Status::DetectorHit) {
	  auto reconstructedPhoton =
	    PhotonReconstructor::ReconstructPhoton(particleTrack,
						   particleTrack.GetPhotonHits().back(),
						   Photon::Radiator::Gas);
	  x += reconstructedPhoton.m_CherenkovAngle;
	  x2 += reconstructedPhoton.m_CherenkovAngle*reconstructedPhoton.m_CherenkovAngle;
	  NumberPhotons++;
	}
      }
    }
    Resolution = NewResolution;
    NewResolution = TMath::Sqrt(x2/NumberPhotons - (x/NumberPhotons)*(x/NumberPhotons));
    NewResolution /= TMath::Sqrt(NumberPhotons);
    if(TMath::Abs(NewResolution - Resolution)/NewResolution < Precision) {
      // Round to nearest "Precision"
      const double AbsolutePrecision = NewResolution*Precision;
      Resolution = std::round(NewResolution/AbsolutePrecision)*AbsolutePrecision;
      return Resolution;
    }
  }
}
