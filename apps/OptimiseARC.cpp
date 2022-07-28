// Martin Duy Tat 1st May 2022
/**
 * OptimiseARC is an application for optimising an ARC cell
 * Run with RunARC ColumnNumber RowNumber <settings name> <settings filename> ...
 * There can be an arbitrary number of settings files added
 */

#include<iostream>
#include<string>
#include<utility>
#include<vector>
#include<fstream>
#include"TRandom.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Settings.h"
#include"RadiatorCell.h"
#include"TrackingVolume.h"
#include"ParticleTrack.h"
#include"ResolutionUtilities.h"

using Vector = ROOT::Math::XYZVector;
using TracksPhotons = std::vector<std::pair<ParticleTrack, std::vector<Photon>>>;

Vector VectorFromSpherical(double R, double CosTheta, double Phi);

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
  std::cout << "Generating tracks and photons...\n";
  const int ParticleID = Settings::GetInt("Particle/ID");;
  const TrackingVolume InnerTracker;
  const int CellsPerRow = 2*Settings::GetDouble("ARCGeometry/CellsPerRow") - 1;
  const double HexagonSize = Settings::GetDouble("ARCGeometry/Length")/CellsPerRow;
  RadiatorCell radiatorCell(std::stoi(std::string(argv[1])),
			    std::stoi(std::string(argv[2])), 
			    HexagonSize);
  const Vector CellPosition = radiatorCell.GetRadiatorPosition();
  const double Radius = Settings::GetDouble("ARCGeometry/Radius");
  const double MomentumMag = Settings::GetDouble("Particle/Momentum");
  const double z_min = CellPosition.X() - HexagonSize/2.0;
  const double z_max = CellPosition.X() + HexagonSize/2.0;
  TracksPhotons ParticlesPhotons;
  const double TotalNumberTracks = Settings::GetInt("General/NumberTracks");
  ParticlesPhotons.reserve(TotalNumberTracks);
  int NumberTracks = 0;
  while(NumberTracks < TotalNumberTracks) {
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
    ParticlesPhotons.push_back(std::make_pair(particleTrack, std::move(Photons)));
    NumberTracks++;
  }
  if(Settings::GetBool("Optimisation/DoFit")) {
    std::cout << "Differential evolution ready, sending off agents...\n";
    ResolutionUtilities::DoFit(radiatorCell, ParticlesPhotons);
    std::cout << "ARC is optimised!\n";
  }
  if(Settings::GetBool("Optimisation/PlotProjections")) {
    std::cout << "Plotting...\n";
    ResolutionUtilities::PlotProjections(radiatorCell, ParticlesPhotons);
    std::cout << "Resolution projections plotted\n";
  }
  return 0;
}

Vector VectorFromSpherical(double R, double CosTheta, double Phi) {
  const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
  const double CosPhi = TMath::Cos(Phi);
  const double SinPhi = TMath::Sin(Phi);
  return Vector{R*CosPhi*SinTheta, R*SinPhi*SinTheta, R*CosTheta};
}


