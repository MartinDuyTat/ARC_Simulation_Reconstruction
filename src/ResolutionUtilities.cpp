// Martin Duy Tat 28th July 2022

#include<omp.h>
#include<fstream>
#include<iostream>
#include<sstream>
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
      bool HitTopWall = false;
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
	} else if(Photon.m_Status == Photon::Status::WallMiss ||
		  Photon.m_Status == Photon::Status::Backwards ||
		  Photon.m_Status == Photon::Status::DetectorMiss) {
	  resolutionStruct.HitTopWall = true;
	}
      }
      return resolutionStruct;
    };
    // Put it all together
    ResolutionStruct Total{};
    std::size_t TracksWithPhotonsHittingWall = 0;
    #pragma omp parallel for num_threads(8)
    for(std::size_t i = 0; i < ParticlesPhotons.size(); i++) {
    ResolutionStruct resolutionStruct = TrackPhotons(ParticlesPhotons[i]);
    #pragma omp critical (Update)
    {
      if(resolutionStruct.HitTopWall) {
	TracksWithPhotonsHittingWall++;
      }
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
      const auto Particles = ParticlesPhotons.size();
      const double Penalty = (1.0*TracksWithPhotonsHittingWall)/Particles;
      return Resolution/TMath::Sqrt(MeanNumberPhotons) + Penalty;
    } 
  }

  double fcn(double MirrorCurvature,
	     double MirrorXPosition,
	     double MirrorZPosition,
	     double DetectorPosition,
	     double DetectorTilt,
	     RadiatorCell &radiatorCell,
	     const TracksPhotons &ParticlesPhotons) {
    radiatorCell.SetMirrorCurvature(MirrorCurvature);
    radiatorCell.SetMirrorXPosition(MirrorXPosition);
    radiatorCell.SetMirrorZPosition(MirrorZPosition);
    radiatorCell.SetDetectorPosition(DetectorPosition);
    radiatorCell.SetDetectorTilt(DetectorTilt);
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
      std::string dummy;
      double Value;
      std::stringstream ss(Line);
      ss >> dummy >> Value;
      Result.push_back(Value);
    }
    File.close();
    auto MinimiseFunctionMirrorCurvature = [&] (double *x, double*) {
      return fcn(x[0], Result[1], Result[2], Result[3], Result[4],
	         radiatorCell, ParticlesPhotons);
    };
    auto MinimiseFunctionXPosition = [&] (double *x, double*) {
      return fcn(Result[0], x[0], Result[2], Result[3], Result[4],
	         radiatorCell, ParticlesPhotons);
    };
    auto MinimiseFunctionZPosition = [&] (double *x, double*) {
      return fcn(Result[0], Result[1], x[0], Result[3], Result[4],
	         radiatorCell, ParticlesPhotons);
    };
    auto MinimiseFunctionDetPosition = [&] (double *x, double*) {
      return fcn(Result[0], Result[1], Result[2], x[0], Result[4],
	         radiatorCell, ParticlesPhotons);
    };
    auto MinimiseFunctionDetTilt = [&] (double *x, double*) {
      return fcn(Result[0], Result[1], Result[2], Result[3], x[0],
	         radiatorCell, ParticlesPhotons);
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

    std::string Name4("Optimisation/DetectorPositionPlot_");
    double DetPosition_min = Settings::GetDouble(Name4 + "min");
    double DetPosition_max = Settings::GetDouble(Name4 + "max");
    TF1 f4("DetectorPosition", MinimiseFunctionDetPosition,
	   DetPosition_min, DetPosition_max, 0);
    std::string Name5("Optimisation/DetectorTiltPlot_");
    double DetTilt_min = Settings::GetDouble(Name5 + "min");
    double DetTilt_max = Settings::GetDouble(Name5 + "max");
    TF1 f5("DetectorTilt", MinimiseFunctionDetTilt,
	   DetTilt_min, DetTilt_max, 0);
    TCanvas c1("c1", "", 1200, 900);
    f1.Draw();
    c1.SaveAs("MirrorCurvatureOptimisation.pdf");
    TCanvas c2("c2", "", 1200, 900);
    f2.Draw();
    c2.SaveAs("XPositionOptimisation.pdf");
    TCanvas c3("c3", "", 1200, 900);
    f3.Draw();
    c3.SaveAs("ZPositionOptimisation.pdf");
    TCanvas c4("c4", "", 1200, 900);
    f4.Draw();
    c4.SaveAs("DetectorPositionOptimisation.pdf");
    TCanvas c5("c4", "", 1200, 900);
    f5.Draw();
    c5.SaveAs("DetectorTiltOptimisation.pdf");
  }

  void DoFit(RadiatorCell &radiatorCell,
	     const TracksPhotons &ParticlesPhotons,
	     int Column,
	     int Row) {
    ResolutionOptimizable resolutionOptimisable(radiatorCell, ParticlesPhotons);
    const int NumberAgents = Settings::GetInt("Optimisation/NumberAgents");
    de::DifferentialEvolution de(resolutionOptimisable, NumberAgents);
    const int Iterations = Settings::GetInt("Optimisation/Iterations");
    de.Optimize(Iterations, true);
    auto Result = de.GetBestAgent();
    auto FixedParameters = resolutionOptimisable.GetFixedParameters();
    std::size_t TotalParameters = Result.size() + FixedParameters.size();
    std::vector<double> AllParameters;
    constexpr std::array<std::string_view, 5> ParameterNames{
      "Curvature",
      "XPosition",
      "ZPosition",
      "DetPosition",
      "DetTilt"};
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
    const double Resolution = fcn(AllParameters[0],
	                          AllParameters[1],
	                          AllParameters[2],
	                          AllParameters[3],
	                          AllParameters[4],
	                          radiatorCell,
	                          ParticlesPhotons);
    std::string ResultFilename = Settings::GetString("Optimisation/Filename");
    std::ofstream File(ResultFilename);
    for(std::size_t i = 0; i < TotalParameters; i++) {
      File << "Radiator_c" << Column << "_r" << Row << "_" << ParameterNames[i] << " ";
      File << AllParameters[i] << "\n";
    }
    File << "\n" << "OptimalResolution: " << Resolution << "\n";
    File.close();
  }

}
