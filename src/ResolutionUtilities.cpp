// Martin Duy Tat 28th July 2022

#include<omp.h>
#include<fstream>
#include<iostream>
#include"TCanvas.h"
#include"TF1.h"
#include"ResolutionUtilities.h"
#include"RadiatorCell.h"
#include"ParticleTrack.h"
#include"Photon.h"
#include"PhotonMapper.h"
#include"PhotonReconstructor.h"
#include"SiPM.h"
#include"DifferentialEvolution.h"
#include"ResolutionOptimizable.h"
#include"Settings.h"

namespace ResolutionUtilities {

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
    const double Resolution =
      TMath::Sqrt(Total.x2/Total.N - (Total.x/Total.N)*(Total.x/Total.N));
    const double MeanNumberPhotons =
      static_cast<double>(Total.N)/ParticlesPhotons.size();
    if(MeanNumberPhotons <= 3.0) {
      return 1000.0;
    } else {
      return Resolution/TMath::Sqrt(MeanNumberPhotons);
    } 
  }

  double fcn(double MirrorCurvature,
	     double MirrorXPosition,
	     double MirrorZPosition,
	     RadiatorCell &radiatorCell,
	     const TracksPhotons &ParticlesPhotons) {
    radiatorCell.SetMirrorCurvature(MirrorCurvature);
    radiatorCell.SetMirrorXPosition(MirrorXPosition);
    radiatorCell.SetMirrorZPosition(MirrorZPosition);
    const double Resolution = CalculateResolution(radiatorCell, ParticlesPhotons);
    return Resolution;
  }

  void PlotProjections(RadiatorCell &radiatorCell,
                       const TracksPhotons &ParticlesPhotons) {
    std::string ResultFilename = Settings::GetString("Optimisation/Filename");
    std::ifstream File(ResultFilename);
    std::vector<double> Result;
    std::string Line;
    while(std::getline(File, Line)) {
      Result.push_back(std::stod(Line));
    }
    File.close();
    auto MinimiseFunctionMirrorCurvature = [&] (double *x, double*) {
      return fcn(x[0], Result[1], Result[2], radiatorCell, ParticlesPhotons);
    };
    auto MinimiseFunctionXPosition = [&] (double *x, double*) {
      return fcn(Result[0], x[0], Result[2], radiatorCell, ParticlesPhotons);
    };
    auto MinimiseFunctionZPosition = [&] (double *x, double*) {
      return fcn(Result[0], Result[1], x[0], radiatorCell, ParticlesPhotons);
    };
    std::string Name1("Optimisation/MirrorCurvaturePlot_");
    double Curvature_min = Settings::GetDouble(Name1 + "min");
    double Curvature_max = Settings::GetDouble(Name1 + "max");
    TF1 f1("MirrorCurvature", MinimiseFunctionMirrorCurvature,
	   Curvature_min, Curvature_max, 0);
    std::string Name2("Optimisation/MirrorXPositionPlot_");
    double XPosition_min = Settings::GetDouble(Name2 + "min");
    double XPosition_max = Settings::GetDouble(Name2 + "max");
    TF1 f2("XPosition", MinimiseFunctionXPosition,
	   XPosition_min, XPosition_max, 0);
    std::string Name3("Optimisation/MirrorZPositionPlot_");
    double ZPosition_min = Settings::GetDouble(Name3 + "min");
    double ZPosition_max = Settings::GetDouble(Name3 + "max");
    TF1 f3("ZPosition", MinimiseFunctionZPosition,
	   ZPosition_min, ZPosition_max, 0);
    TCanvas c1("c1", "", 1200, 900);
    f1.Draw();
    c1.SaveAs("MirrorCurvatureOptimisation.pdf");
    TCanvas c2("c2", "", 1200, 900);
    f2.Draw();
    c2.SaveAs("XPositionOptimisation.pdf");
    TCanvas c3("c3", "", 1200, 900);
    f3.Draw();
    c3.SaveAs("ZPositionOptimisation.pdf");
  }

  void DoFit(RadiatorCell &radiatorCell,
	     const TracksPhotons &ParticlesPhotons) {
    ResolutionOptimizable resolutionOptimisable(radiatorCell, ParticlesPhotons);
    const int NumberAgents = Settings::GetInt("Optimisation/NumberAgents");
    de::DifferentialEvolution de(resolutionOptimisable, NumberAgents);
    const int Iterations = Settings::GetInt("Optimisation/Iterations");
    de.Optimize(Iterations, true);
    auto Result = de.GetBestAgent();
    auto FixedParameters = resolutionOptimisable.GetFixedParameters();
    std::size_t TotalParameters = Result.size() + FixedParameters.size();
    std::vector<double> AllParameters;
    std::size_t j = 0;
    for(std::size_t i = 0; i < TotalParameters; i++) {
      auto iter = FixedParameters.find(i);
      if(iter != FixedParameters.end()) {
	AllParameters.push_back(iter->second);
      } else {
	AllParameters.push_back(Result[j]);
	j++;
      }
    }
    std::string ResultFilename = Settings::GetString("Optimisation/Filename");
    std::ofstream File(ResultFilename);
    for(const auto &Parameter : AllParameters) {
      File << Parameter << "\n";
    }
    File.close();
  }

}
