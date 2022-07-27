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
#include<numeric>
#include"TRandom.h"
#include"TF1.h"
#include"TCanvas.h"
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
#include"SiPM.h"

using Vector = ROOT::Math::XYZVector;
using TracksPhotons = std::vector<std::pair<ParticleTrack, std::vector<Photon>>>;

Vector VectorFromSpherical(double R, double CosTheta, double Phi);

double CalculateResolution(const RadiatorCell &radiatorCell,
			   const TracksPhotons &ParticlesPhotons);
			   

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
  const double HexagonSize = Settings::GetDouble("ARCGeometry/Length")/(2*CellsPerRow);
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
  auto MinimiseFunction = [&] (const double *x) {
    radiatorCell.SetMirrorCurvature(x[0]);
    radiatorCell.SetMirrorXPosition(x[1]);
    radiatorCell.SetMirrorZPosition(x[2]);
    const double Resolution = CalculateResolution(radiatorCell, ParticlesPhotons);
    return Resolution;
  };
  ROOT::Math::Functor fcn(MinimiseFunction, 3);
  ROOT::Minuit2::Minuit2Minimizer Minimiser(ROOT::Minuit2::kSimplex);
  Minimiser.SetPrintLevel(4);
  Minimiser.SetPrecision(0.0001);
  Minimiser.SetFunction(fcn);
  Minimiser.SetVariable(0, "MirrorCurvature", 0.37, 0.01);
  Minimiser.SetVariableLimits(0, 0.33, 0.40);
  Minimiser.SetVariable(1, "Mirror_xPosition", 0.01, 0.005);
  Minimiser.SetVariableLimits(1, -0.1, 0.1);
  //Minimiser.FixVariable(1);
  Minimiser.SetVariable(2, "Mirror_zPosition", 0.0, 0.01);
  Minimiser.SetVariableLimits(2, -0.05, 0.05);
  std::cout << "Starting minimisation...\n";
  Minimiser.Minimize();
  std::cout << "Status: " << Minimiser.Status() << "\n";
  std::cout << "Covariance matrix status: " << Minimiser.CovMatrixStatus() << "\n";
  std::cout << "ARC is optimised!\n";
  if(Settings::GetBool("General/PlotProjections")) {
    std::cout << "Plotting...\n";
    auto Result = Minimiser.X();
    auto MinimiseFunctionMirrorCurvature = [&] (double *x, double*) {
      double xx[3] = {x[0], Result[1], Result[2]};
      return MinimiseFunction(xx);
    };
    auto MinimiseFunctionXPosition = [&] (double *x, double*) {
      double xx[3] = {Result[0], x[0], Result[2]};
      return MinimiseFunction(xx);
    };
    auto MinimiseFunctionZPosition = [&] (double *x, double*) {
      double xx[3] = {Result[0], Result[1], x[0]};
      return MinimiseFunction(xx);
    };
    TF1 f1("MirrorCurvature", MinimiseFunctionMirrorCurvature, 0.34, 0.375, 0);
    TF1 f2("XPosition", MinimiseFunctionXPosition, -0.05, 0.05, 0);
    TF1 f3("ZPosition", MinimiseFunctionZPosition, -0.002, 0.0045, 0);
    TCanvas c1("c1", "", 1200, 900);
    f1.Draw();
    c1.SaveAs("MirrorCurvatureOptimisation.pdf");
    TCanvas c2("c2", "", 1200, 900);
    f2.Draw();
    c2.SaveAs("XPositionOptimisation.pdf");
    TCanvas c3("c3", "", 1200, 900);
    f3.Draw();
    c3.SaveAs("ZPositionOptimisation.pdf");
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

double CalculateResolution(const RadiatorCell &radiatorCell,
			   const TracksPhotons &ParticlesPhotons) {
  // Struct containing sum of Cherenkov angle, sum of square of Cherenkov
  // angle and total number of photons
  struct ResolutionStruct {
    // x is the Cherenkov angle, x2 is the Cherenkov angle squared
    double x = 0.0;
    double x2 = 0.0;
    int N = 0;
  };
  // Lambda for tracking photons and calculating the Cherenkov angles
  auto TrackPhotons = [&radiatorCell] (const auto &a) {
    const auto particleTrack = a.first;
    const auto Photons = a.second;
    ResolutionStruct resolutionStruct;
    for(auto Photon : Photons) {
      auto photonHit = PhotonMapper::TracePhoton(Photon);
      if(Photon.m_Status == Photon::Status::DetectorHit) {
	auto reconstructedPhoton =
	  PhotonReconstructor::ReconstructPhoton(particleTrack,
						 *photonHit,
						 Photon::Radiator::Gas);
	const double CherenkovAngle = TMath::ACos(reconstructedPhoton.m_CosCherenkovAngle);
	resolutionStruct.x += CherenkovAngle;
	resolutionStruct.x2 += CherenkovAngle*CherenkovAngle;
	resolutionStruct.N++;
      }
    }
    return resolutionStruct;
  };
  // Put it all together
  ResolutionStruct Total{};
  #pragma omp parallel for num_threads(8)
  for(std::size_t i = 0; i < ParticlesPhotons.size(); i++) {
    ResolutionStruct resolutionStruct = TrackPhotons(ParticlesPhotons[i]);
    #pragma omp critical (Update)
    {
    Total.x += resolutionStruct.x;
    Total.x2 += resolutionStruct.x2;
    Total.N += resolutionStruct.N;
    }
  }
  const double Resolution = TMath::Sqrt(Total.x2/Total.N
			  - (Total.x/Total.N)*(Total.x/Total.N));
  const double MeanNumberPhotons = static_cast<double>(Total.N)
                                   /ParticlesPhotons.size();
  if(MeanNumberPhotons <= 0.0) {
    return 1000.0;
  } else {
    return Resolution/TMath::Sqrt(MeanNumberPhotons);
  }
}
