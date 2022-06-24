// Martin Duy Tat 1st May 2022
/**
 * RunARC is an application for running ARC simulations and reconstructions
 * Run with RunARC <run mode> <settings name> <settings filename> ...
 * There can be an arbitrary number of settings files added
 * Possible run modes:
 * "SingleTrack" Send a single charged track and generate a photon hit map
 * "CherenkovAngleResolution" Generate photons from a single track, reconstruct the photons and study the Cherenkov angle resolution
 */

#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include"TFile.h"
#include"TTree.h"
#include"TMath.h"
#include"TRandom.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"ParticleTrack.h"
#include"TrackingVolume.h"
#include"Photon.h"
#include"PhotonMapper.h"
#include"PhotonReconstructor.h"
#include"RadiatorArray.h"
#include"Settings.h"
#include"EventDisplay.h"

using Vector = ROOT::Math::XYZVector;

Vector VectorFromSpherical(double R, double Theta, double Phi);

int main(int argc, char *argv[]) {
  if(argc%2 != 0) {
    return 0;
  }
  std::cout << "Welcome to the ARC simulation and reconstruction\n";
  for(int i = 2; i < argc; i += 2) {
    const std::string SettingsName = argv[i];
    const std::string SettingsFilename = argv[i + 1];
    Settings::AddSettings(SettingsName, SettingsFilename);
    std::cout << "Added " << SettingsName << " settings from " << SettingsFilename << "\n";
  }
  const std::string RunMode(argv[1]);
  if(Settings::Exists("General/Seed")) {
    gRandom->SetSeed(Settings::GetInt("General/Seed"));
  }
  EventDisplay eventDisplay;
  const TrackingVolume InnerTracker;
  eventDisplay.AddObject(InnerTracker.DrawARCGeometry());
  RadiatorArray radiatorArray;
  eventDisplay.AddObject(radiatorArray.DrawRadiatorArray());
  if(RunMode == "SingleTrack") {
    std::cout << "Run mode: Single track\n";
    const Vector Momentum = VectorFromSpherical(Settings::GetDouble("Particle/Momentum"),
						Settings::GetDouble("Particle/Theta"),
						Settings::GetDouble("Particle/Phi"));
    const int ParticleID = Settings::GetInt("Particle/ID");;
    ParticleTrack particleTrack(Momentum, ParticleID);
    particleTrack.TrackThroughTracker(InnerTracker);
    particleTrack.ConvertToRadiatorCoordinates(radiatorArray);
    particleTrack.TrackThroughRadiatorCell();
    auto PhotonsAerogel = particleTrack.GeneratePhotonsFromAerogel();
    auto PhotonsGas = particleTrack.GeneratePhotonsFromGas();
    for(auto &photon : PhotonsAerogel) {
      PhotonMapper::TracePhoton(photon);
    }
    for(auto &photon : PhotonsGas) {
      PhotonMapper::TracePhoton(photon);
    }
    radiatorArray[0].m_Detector.PlotHits("PhotonHits.png");
  } else if(RunMode == "CherenkovAngleResolution") {
    std::cout << "Run mode: Cherenkov angle resolution\n";
    TFile CherenkovFile("CherenkovFile.root", "RECREATE");
    TTree CherenkovTree("CherenkovTree", "");
    double CherenkovAngle_Reco_TrueEmissionPoint, CherenkovAngle_Reco, CherenkovAngle_True, PhotonEnergy;
    int DetectorHit, RadiatorNumber;
    CherenkovTree.Branch("CherenkovAngle_Reco_TrueEmissionPoint", &CherenkovAngle_Reco_TrueEmissionPoint);
    CherenkovTree.Branch("CherenkovAngle_Reco", &CherenkovAngle_Reco);
    CherenkovTree.Branch("CherenkovAngle_True", &CherenkovAngle_True);
    CherenkovTree.Branch("PhotonEnergy", &PhotonEnergy);
    CherenkovTree.Branch("DetectorHit", &DetectorHit);
    CherenkovTree.Branch("RadiatorNumber", &RadiatorNumber);
    const int NumberTracks = Settings::GetInt("General/NumberTracks");
    const std::vector<int> TracksToDraw = Settings::GetIntVector("General/TrackToDraw");
    const bool DrawThetaMiss = Settings::GetBool("General/DrawThetaMiss");
    const bool DrawPhiMiss = Settings::GetBool("General/DrawPhiMiss");
    for(int i = 0; i < NumberTracks; i++) {
      const bool DrawThisTrack = std::find(TracksToDraw.begin(), TracksToDraw.end(), i) != TracksToDraw.end();
      const double Phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
      const double Theta = gRandom->Uniform(Settings::GetDouble("Particle/Theta_min"), Settings::GetDouble("Particle/Theta_max"));;
      const Vector Momentum = VectorFromSpherical(Settings::GetDouble("Particle/Momentum"), Theta, Phi);
      const int ParticleID = Settings::GetInt("Particle/ID");;
      ParticleTrack particleTrack(Momentum, ParticleID);
      particleTrack.TrackThroughTracker(InnerTracker);
      particleTrack.ConvertToRadiatorCoordinates(radiatorArray);
      particleTrack.TrackThroughRadiatorCell();
      if(DrawThisTrack) {
	eventDisplay.AddObject(particleTrack.DrawParticleTrack());
      }
      auto Photons = particleTrack.GeneratePhotonsFromGas();
      for(auto &Photon : Photons) {
	CherenkovAngle_True = Photon.m_CherenkovAngle;
	PhotonMapper::TracePhoton(Photon);
	if(DrawThisTrack) {
	  eventDisplay.AddObject(Photon.DrawPhotonPath());
	}
	if(!Photon.m_MirrorHitPosition) {
	  if(DrawThetaMiss && Photon.m_Status == Photon::Status::MissedTheta) {
	    eventDisplay.AddObject(Photon.DrawPhotonPath());
	  } else if(DrawPhiMiss && Photon.m_Status == Photon::Status::MissedPhi) {
	    eventDisplay.AddObject(Photon.DrawPhotonPath());
	  }
	  continue;
	}
	auto reconstructedPhoton = PhotonReconstructor::ReconstructPhoton(particleTrack, particleTrack.GetPhotonHits().back(), Photon::Radiator::Gas);
	PhotonEnergy = Photon.m_Energy;
	DetectorHit = Photon.m_Status == Photon::Status::DetectorHit ? 1 : 0;
	RadiatorNumber = Photon.m_RadiatorCell->GetCellNumber();
	CherenkovAngle_Reco_TrueEmissionPoint = reconstructedPhoton.m_CherenkovAngle_TrueEmissionPoint;
	CherenkovAngle_Reco = reconstructedPhoton.m_CherenkovAngle;
	CherenkovTree.Fill();
      }
    }
    eventDisplay.DrawEventDisplay("EventDisplay.pdf");
    CherenkovTree.Write();
    CherenkovFile.Close();
  }
  return 0;
}

Vector VectorFromSpherical(double R, double Theta, double Phi) {
  const double CosTheta = TMath::Cos(Theta);
  const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
  const double CosPhi = TMath::Cos(Phi);
  const double SinPhi = TMath::Sqrt(1.0 - CosPhi*CosPhi);
  return Vector{R*CosPhi*SinTheta, R*SinPhi*SinTheta, R*CosTheta};
}
