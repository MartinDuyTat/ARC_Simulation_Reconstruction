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
#include<array>
#include<memory>
#include"TRandom.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Settings.h"
#include"RadiatorArray.h"
#include"BarrelRadiatorArray.h"
#include"EndCapRadiatorArray.h"
#include"RadiatorCell.h"
#include"HalfRadiatorCell.h"
#include"TrackingVolume.h"
#include"ParticleTrack.h"
#include"ResolutionUtilities.h"
#include"Utilities.h"

using Vector = ROOT::Math::XYZVector;
using Tracks = std::vector<ParticleTrack>;

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
  const std::size_t Seed = Settings::GetSizeT("General/Seed");
  gRandom->SetSeed(Seed);
  std::cout << "Generating tracks...\n";
  const int ParticleID = Settings::GetInt("Particle/ID");;
  const TrackingVolume InnerTracker;
  const std::size_t Column = static_cast<std::size_t>(std::stoi(std::string(argv[1])));
  const std::size_t Row = static_cast<std::size_t>(std::stoi(std::string(argv[2])));
  std::cout << "Tracks ready\n";
  std::unique_ptr<RadiatorArray> radiatorArray;
  const std::string BarrelOrEndcap = Settings::GetString("General/BarrelOrEndcap");
  if(BarrelOrEndcap == "Barrel") {
    radiatorArray = std::make_unique<BarrelRadiatorArray>();
  } else if(BarrelOrEndcap == "EndCap") {
    radiatorArray = std::make_unique<EndCapRadiatorArray>();
  } else {
    return 0;
  }
  RadiatorCell *radiatorCell = radiatorArray->GetRadiatorCell(Column, Row);
  const double z_min = Settings::GetDouble("Particle/z_min");
  const double z_max = Settings::GetDouble("Particle/z_max");
  Tracks Particles;
  const std::size_t TotalNumberTracks = Settings::GetSizeT("General/NumberTracks");
  Particles.reserve(TotalNumberTracks);
  std::size_t NumberTracks = 0;
  while(NumberTracks < TotalNumberTracks) {
    auto GetMomentum = [=] () {
      if(BarrelOrEndcap == "Barrel") {
	return Utilities::GenerateRandomBarrelTrackZRange(z_min, z_max);
      } else if(BarrelOrEndcap == "EndCap") {
	return Utilities::GenerateRandomEndCapTrack();
      } else {
	return Vector(0.0, 0.0, 0.0);
      }
    };
    const Vector Momentum = GetMomentum();
    ParticleTrack particleTrack(ParticleID, Momentum, NumberTracks);
    particleTrack.TrackThroughTracker(InnerTracker);
    Particles.push_back(particleTrack);
    NumberTracks++;
  }
  if(Settings::GetBool("Optimisation/SinglePoints")) {
    while(true) {
      std::cout << "Input parameters:\n";
      std::array<double, 5> x;
      for(double &xx : x) {
	std::cin >> xx;
      }
      double Value = ResolutionUtilities::fcn(x[0], x[1], x[2], x[3], x[4],
					      *radiatorCell, *radiatorArray,
					      Particles, Seed, false);
      std::cout << "fcn = " << Value << "\n";
      std::cout << "Calculate new value?(y/n)\n";
      std::string Answer;
      std::cin >> Answer;
      if(Answer != "y") {
	return 0;
      }
    }
  }
  if(Settings::GetBool("Optimisation/DoFit")) {
    std::cout << "Differential evolution ready, sending off agents...\n";
    ResolutionUtilities::DoFit(*radiatorCell, *radiatorArray, Particles, Column, Row);
    std::cout << "ARC is optimised!\n";
  }
  if(Settings::GetBool("Optimisation/PlotProjections")) {
    std::cout << "Plotting...\n";
    ResolutionUtilities::PlotProjections(*radiatorCell, *radiatorArray, Particles);
    std::cout << "Resolution projections plotted\n";
  }
  return 0;
}
