// Martin Duy Tat 28th July 2022

#include<omp.h>
#include<fstream>
#include<iostream>
#include<sstream>
#include"TCanvas.h"
#include"TF1.h"
#include"TLine.h"
#include"TRandom.h"
#include"TPad.h"
#include"TH1.h"
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

  double CalculateResolution(const Tracks &Particles,
			     bool IncludeCentrePenalty) {
    // Struct containing sum of Cherenkov angle, sum of square of Cherenkov
    // angle and total number of photons
    struct ResolutionStruct {
      // x is the Cherenkov angle, x2 is the Cherenkov angle squared
      // CentreHitDistance is the distance between photon hit and detector centre
      double x = 0.0;
      double x2 = 0.0;
      int N = 0;
      bool HitTopWall = false;
      Vector CentreHitDistance{0.0, 0.0, 0.0};
    };
    // Lambda for tracking photons and calculating the Cherenkov angles
    auto TrackPhotons = [] (auto particleTrack) {
      particleTrack.TrackThroughRadiatorCell();
      auto Photons = particleTrack.GetParticleLocation() != ParticleTrack::Location::Mirror ?
                     std::vector<Photon>() : particleTrack.GeneratePhotonsFromGas();
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
	  resolutionStruct.CentreHitDistance += photonHit->m_CentreHitDistance;
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
    std::size_t PhotonsHitWallTracks = 0;
    #pragma omp parallel for num_threads(8)
    for(std::size_t i = 0; i < Particles.size(); i++) {
    ResolutionStruct resolutionStruct = TrackPhotons(Particles[i]);
    #pragma omp critical (Update)
    {
      if(resolutionStruct.HitTopWall) {
	PhotonsHitWallTracks++;
      }
      Total.x += resolutionStruct.x;
      Total.x2 += resolutionStruct.x2;
      Total.N += resolutionStruct.N;
      Total.CentreHitDistance += resolutionStruct.CentreHitDistance;
    }
    }
    // 0.0003 is the pixel size angular resolution
    const double Resolution =
      TMath::Sqrt(Total.x2/Total.N - (Total.x/Total.N)*(Total.x/Total.N)
	        + 0.0003*0.0003);
    const double MeanNumberPhotons =
      static_cast<double>(Total.N)/static_cast<double>(Particles.size());
    if(MeanNumberPhotons <= 1.0) {
      return 1000.0;
    } else {
      const auto NumberParticles = Particles.size();
      auto GetFinalResolution = [&] (bool IncludeCentrePenalty) {
        // Penalty when hitting the upper wall
        const double WallPenalty = (10.0*static_cast<double>(PhotonsHitWallTracks))/
	                           static_cast<double>(NumberParticles);
        const double ResolutionWithPenalty = Resolution/TMath::Sqrt(MeanNumberPhotons)
	                                   + WallPenalty;
	if(IncludeCentrePenalty) {
          // Penalty when far from the mirror centre
          const Vector AverageRingPosition = Total.CentreHitDistance*(1.0/Total.N);
          const double CentrePenalty = 0.01*TMath::Sqrt(AverageRingPosition.Mag2());
          return ResolutionWithPenalty + CentrePenalty;
        } else {
          return ResolutionWithPenalty;
        }
      };
      const double FinalResolution = GetFinalResolution(IncludeCentrePenalty);
      return FinalResolution;
    } 
  }

  double fcn(double MirrorCurvature,
	     double MirrorXPosition,
	     double MirrorZPosition,
	     double DetectorPosition,
	     double DetectorTilt,
	     RadiatorCell &radiatorCell,
	     const Tracks &Particles,
	     std::size_t Seed,
	     bool IncludeCentrePenalty) {
    radiatorCell.SetMirrorCurvature(MirrorCurvature);
    radiatorCell.SetMirrorXPosition(MirrorXPosition);
    radiatorCell.SetMirrorZPosition(MirrorZPosition);
    radiatorCell.SetDetectorPosition(DetectorPosition);
    radiatorCell.SetDetectorTilt(DetectorTilt);
    if(!radiatorCell.IsDetectorInsideCell()) {
      return 1000.0;
    }
    gRandom->SetSeed(Seed);
    const double Resolution = CalculateResolution(Particles,
	                                          IncludeCentrePenalty);
    return Resolution;
  }

  void PlotProjections(RadiatorCell &radiatorCell,
                       const Tracks &Particles) {
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
    const std::size_t Seed = Settings::GetSizeT("General/Seed");
    auto MinimiseFunctionMirrorCurvature = [&] (double *x, double*) {
      return fcn(x[0], Result[1], Result[2], Result[3], Result[4],
	         radiatorCell, Particles, Seed, false);
    };
    auto MinimiseFunctionXPosition = [&] (double *x, double*) {
      return fcn(Result[0], x[0], Result[2], Result[3], Result[4],
	         radiatorCell, Particles, Seed, false);
    };
    auto MinimiseFunctionZPosition = [&] (double *x, double*) {
      return fcn(Result[0], Result[1], x[0], Result[3], Result[4],
	         radiatorCell, Particles, Seed, false);
    };
    auto MinimiseFunctionDetPosition = [&] (double *x, double*) {
      return fcn(Result[0], Result[1], Result[2], x[0], Result[4],
	         radiatorCell, Particles, Seed, false);
    };
    auto MinimiseFunctionDetTilt = [&] (double *x, double*) {
      return fcn(Result[0], Result[1], Result[2], Result[3], x[0],
	         radiatorCell, Particles, Seed, false);
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
    auto GetSolutionLine = [] (double Result, double ymax = 0.0) {
      TLine Solution(Result, 0.0, Result, 0.9*ymax);
      Solution.SetLineWidth(3);
      return Solution;
    };
    TCanvas c1("c1", "", 1200, 900);
    f1.SetTitle("Mirror curvature;Curvature (m);Resolution (rad)");
    f1.Draw();
    auto Solution1 = GetSolutionLine(Result[0], f1.GetHistogram()->GetMaximum());
    Solution1.Draw("SAME");
    c1.SaveAs("MirrorCurvatureOptimisation.pdf");
    TCanvas c2("c2", "", 1200, 900);
    f2.SetTitle("Horizontal mirror position;x (m);Resolution (rad)");
    f2.Draw();
    auto Solution2 = GetSolutionLine(Result[1], f2.GetHistogram()->GetMaximum());
    Solution2.Draw("SAME");
    c2.SaveAs("XPositionOptimisation.pdf");
    TCanvas c3("c3", "", 1200, 900);
    f3.SetTitle("Vertical mirror position;z (m);Resolution (rad)");
    f3.Draw();
    auto Solution3 = GetSolutionLine(Result[2], f3.GetHistogram()->GetMaximum());
    Solution3.Draw("SAME");
    c3.SaveAs("ZPositionOptimisation.pdf");
    TCanvas c4("c4", "", 1200, 900);
    f4.SetTitle("Horizontal detector position;x (m);Resolution (rad)");
    f4.Draw();
    auto Solution4 = GetSolutionLine(Result[3], f4.GetHistogram()->GetMaximum());
    Solution4.Draw("SAME");
    c4.SaveAs("DetectorPositionOptimisation.pdf");
    TCanvas c5("c4", "", 1200, 900);
    f5.SetTitle("Detector tilt angle;#theta (^{o});Resolution (rad)");
    f5.Draw();
    auto Solution5 = GetSolutionLine(Result[4], f5.GetHistogram()->GetMaximum());
    Solution5.Draw("SAME");
    c5.SaveAs("DetectorTiltOptimisation.pdf");
  }

  void DoFit(RadiatorCell &radiatorCell,
	     const Tracks &Particles,
	     int Column,
	     int Row) {
    ResolutionOptimizable resolutionOptimisable(radiatorCell, Particles);
    const std::size_t NumberAgents = Settings::GetSizeT("Optimisation/NumberAgents");
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
    std::string ResultFilename = Settings::GetString("Optimisation/Filename");
    std::ofstream File(ResultFilename);
    const std::string BarrelOrEndcap = Settings::GetString("General/BarrelOrEndcap");
    const std::string Prefix = BarrelOrEndcap == "Barrel" ? "" : "EndCap";
    for(std::size_t i = 0; i < TotalParameters; i++) {
      File << Prefix;
      File << "Radiator_c" << Column << "_r" << Row << "_" << ParameterNames[i] << " ";
      File << AllParameters[i] << "\n";
    }
    //File << "\n" << "OptimalResolution: " << de.GetBestCost()*1000.0 << " mrad" << "\n";
    File.close();
  }

}
