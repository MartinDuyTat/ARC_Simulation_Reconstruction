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
#include<stdexcept>
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
#include"SiPM.h"
#include"BarrelRadiatorArray.h"
#include"EndCapRadiatorArray.h"
#include"Utilities.h"

using Vector = ROOT::Math::XYZVector;

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
  std::unique_ptr<RadiatorArray> radiatorArray;
  const std::string BarrelOrEndcap = Settings::GetString("General/BarrelOrEndcap");
  if(BarrelOrEndcap == "Barrel") {
    radiatorArray = std::make_unique<BarrelRadiatorArray>();
  } else if(BarrelOrEndcap == "EndCap") {
    radiatorArray = std::make_unique<EndCapRadiatorArray>();
  } else {
    return 0;
  } 
  eventDisplay.AddObject(radiatorArray->DrawRadiatorArray());
  if(RunMode == "SingleTrack") {
    std::cout << "Run mode: Single track\n";
    const Vector Momentum = Utilities::VectorFromSpherical(
      Settings::GetDouble("Particle/Momentum"),
      Settings::GetDouble("Particle/CosTheta"),
      Settings::GetDouble("Particle/Phi"));
    const int ParticleID = Settings::GetInt("Particle/ID");;
    ParticleTrack particleTrack(ParticleID, Momentum);
    particleTrack.TrackThroughTracker(InnerTracker);
    particleTrack.FindRadiator(*radiatorArray);
    particleTrack.ConvertToRadiatorCoordinates();
    particleTrack.TrackThroughRadiatorCell();
    auto PhotonsAerogel = particleTrack.GeneratePhotonsFromAerogel();
    auto PhotonsGas = particleTrack.GeneratePhotonsFromGas();
    std::vector<PhotonHit> photonHits;
    for(auto &photon : PhotonsAerogel) {
      auto photonHit = PhotonMapper::TracePhoton(photon);
      if(!photonHit) {
	photonHits.push_back(*photonHit);
      }
    }
    for(auto &photon : PhotonsGas) {
      auto photonHit = PhotonMapper::TracePhoton(photon);
      if(!photonHit) {
	photonHits.push_back(*photonHit);
      }
    }
    (*radiatorArray)(0, 0)->GetDetector().PlotHits("PhotonHits.pdf", photonHits);
  } else if(RunMode == "CherenkovAngleResolution") {
    std::cout << "Run mode: Cherenkov angle resolution\n";
    TFile CherenkovFile("CherenkovFile.root", "RECREATE");
    TTree CherenkovTree("CherenkovTree", "");
    double CherenkovAngle_Reco_TrueEmissionPoint[200], CherenkovAngle_Reco[200],
           CherenkovAngle_True[200], PhotonEnergy[200];
    double Entrance_x, Entrance_y, Entrance_z;
    double CosTheta, Phi;
    int NumberPhotons = 0;
    int RadiatorRowNumber, RadiatorColumnNumber, TrackNumber;
    ParticleTrack::Location ParticleLocation;
    Photon::Status PhotonStatus[200];
    CherenkovTree.Branch("NumberPhotons", &NumberPhotons);
    CherenkovTree.Branch("CherenkovAngle_Reco_TrueEmissionPoint",
			 &CherenkovAngle_Reco_TrueEmissionPoint,
			 "CherenkovAngle_Reco_TrueEmissionPoint[NumberPhotons]/D");
    CherenkovTree.Branch("CherenkovAngle_Reco",
			 &CherenkovAngle_Reco,
			 "CherenkovAngle_Reco[NumberPhotons]/D");
    CherenkovTree.Branch("CherenkovAngle_True",
			 &CherenkovAngle_True,
			 "CherenkovAngle_True[NumberPhotons]/D");
    CherenkovTree.Branch("PhotonEnergy", &PhotonEnergy, "PhotonEnergy[NumberPhotons]/D");
    CherenkovTree.Branch("RadiatorRowNumber",&RadiatorRowNumber);
    CherenkovTree.Branch("RadiatorColumnNumber", &RadiatorColumnNumber);
    CherenkovTree.Branch("TrackNumber", &TrackNumber);
    CherenkovTree.Branch("ParticleLocation", &ParticleLocation, "ParticleLocation/I");
    CherenkovTree.Branch("PhotonStatus", &PhotonStatus, "PhotonStatus[NumberPhotons]/I");
    CherenkovTree.Branch("Entrance_x", &Entrance_x);
    CherenkovTree.Branch("Entrance_y", &Entrance_y);
    CherenkovTree.Branch("Entrance_z", &Entrance_z);
    CherenkovTree.Branch("CosTheta", &CosTheta);
    CherenkovTree.Branch("Phi", &Phi);
    const int NumberTracks = Settings::GetInt("General/NumberTracks");
    const bool DrawAllTracks = Settings::GetBool("General/DrawAllTracks");
    const std::vector<int> TracksToDraw = Settings::GetIntVector("General/TrackToDraw");
    const bool DrawMissPhoton = Settings::GetBool("General/DrawMissPhoton");
    for(int i = 0; i < NumberTracks; i++) {
      NumberPhotons = 0;
      TrackNumber = i;
      RadiatorRowNumber = -1;
      RadiatorColumnNumber = -1;
      auto GetMomentum = [&] () {
	if(BarrelOrEndcap == "Barrel") {
	  return Utilities::GenerateRandomBarrelTrack(CosTheta, Phi);
	} else if(BarrelOrEndcap == "EndCap") {
	  return Utilities::GenerateRandomEndCapTrack();
	} else {
	  return Vector(0.0, 0.0, 0.0);
	}
      };
      const Vector Momentum = GetMomentum();
      const Vector Position(0.0, 0.0, 0.0);
      const int ParticleID = Settings::GetInt("Particle/ID");;
      ParticleTrack particleTrack(ParticleID, Momentum, Position);
      particleTrack.TrackThroughTracker(InnerTracker);
      if(!particleTrack.FindRadiator(*radiatorArray)) {
	continue;
      }
      Phi = particleTrack.GetPosition().Phi();
      RadiatorRowNumber = particleTrack.GetRadiatorRowNumber();
      RadiatorColumnNumber = particleTrack.GetRadiatorColumnNumber();
      auto IsTrackDraw = [&] () {
	if(RadiatorRowNumber != Settings::GetInt("EventDisplay/RowToDraw")) {
	  return false;
	}
	if(DrawAllTracks) {
	  return true;
	}
	const auto iter = std::find(TracksToDraw.begin(), TracksToDraw.end(), i);
	if(iter == TracksToDraw.end()) {
	  return false;
	}
	return true;
      };
      const bool DrawThisTrack = IsTrackDraw();
      particleTrack.ConvertToRadiatorCoordinates();
      auto EntranceWindowPosition = particleTrack.GetEntranceWindowPosition();
      Entrance_x = EntranceWindowPosition.X();
      Entrance_y = EntranceWindowPosition.Y();
      Entrance_z = EntranceWindowPosition.Z();
      if(particleTrack.GetParticleLocation() != ParticleTrack::Location::EntranceWindow) {
	ParticleLocation = particleTrack.GetParticleLocation();
	CherenkovTree.Fill();
	continue;
      }
      particleTrack.TrackThroughRadiatorCell();
      if(particleTrack.GetParticleLocation() != ParticleTrack::Location::Mirror) {
	ParticleLocation = particleTrack.GetParticleLocation();
	CherenkovTree.Fill();
	continue;
      }
      ParticleLocation = particleTrack.GetParticleLocation();
      if(DrawThisTrack) {
	eventDisplay.AddObject(particleTrack.DrawParticleTrack());
      }
      auto Photons = particleTrack.GeneratePhotonsFromGas();
      for(auto &Photon : Photons) {
	CherenkovAngle_True[NumberPhotons] = TMath::ACos(Photon.m_CosCherenkovAngle);
	auto photonHit = PhotonMapper::TracePhoton(Photon);
	if(DrawThisTrack) {
	  eventDisplay.AddObject(Photon.DrawPhotonPath());
	}
	if(!Photon.m_MirrorHitPosition) {
	  if(DrawMissPhoton && Photon.m_Status == Photon::Status::MirrorMiss) {
	    eventDisplay.AddObject(Photon.DrawPhotonPath());
	  }
	  CherenkovAngle_Reco_TrueEmissionPoint[NumberPhotons] = -1.0;
	  CherenkovAngle_Reco[NumberPhotons] = -1.0;
	} else {
	  auto reconstructedPhoton =
	    PhotonReconstructor::ReconstructPhoton(particleTrack,
						   *photonHit,
						   Photon::Radiator::Gas);
	  CherenkovAngle_Reco_TrueEmissionPoint[NumberPhotons] =
	    TMath::ACos(reconstructedPhoton.m_CosCherenkovAngle_TrueEmissionPoint);
	  CherenkovAngle_Reco[NumberPhotons] =
	    TMath::ACos(reconstructedPhoton.m_CosCherenkovAngle);
	}
	PhotonEnergy[NumberPhotons] = Photon.m_Energy;
	PhotonStatus[NumberPhotons] = Photon.m_Status;
	NumberPhotons++;
      }
      CherenkovTree.Fill();
    }
    std::string EventDisplayFilename("EventDisplay");
    if(Settings::GetInt("EventDisplay/RowToDraw") == 1) {
      EventDisplayFilename += "_MainRow";
    } else if(Settings::GetInt("EventDisplay/RowToDraw") == 2) {
      EventDisplayFilename += "_UpperRow";
    }
    EventDisplayFilename += ".pdf";
    eventDisplay.DrawEventDisplay(EventDisplayFilename);
    CherenkovTree.Write();
    CherenkovFile.Close();
  }
  return 0;
}
