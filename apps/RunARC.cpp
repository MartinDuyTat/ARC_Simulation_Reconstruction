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
    gRandom->SetSeed(Settings::GetSizeT("General/Seed"));
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
    particleTrack.TrackThroughAerogel();
    particleTrack.TrackThroughGasToMirror();
    if(particleTrack.GetParticleLocation() != ParticleTrack::Location::Mirror) {
      particleTrack.TrackToNextCell(*radiatorArray);
    }
    auto PhotonsAerogel = particleTrack.GeneratePhotonsFromAerogel();
    auto PhotonsGas = particleTrack.GeneratePhotonsFromGas();
    std::vector<PhotonHit> photonHits;
    for(auto &photon : PhotonsAerogel) {
      auto photonHit = PhotonMapper::TracePhoton(photon, *radiatorArray);
      if(photonHit) {
	photonHits.push_back(*photonHit);
      }
    }
    for(auto &photon : PhotonsGas) {
      auto photonHit = PhotonMapper::TracePhoton(photon, *radiatorArray);
      if(photonHit) {
	photonHits.push_back(*photonHit);
      }
    }
    (*radiatorArray)(0, 0)->GetDetector().PlotHits("PhotonHits.pdf", photonHits);
  } else if(RunMode == "CherenkovAngleResolution") {
    std::cout << "Run mode: Cherenkov angle resolution\n";
    TFile CherenkovFile("CherenkovFile.root", "RECREATE");
    TTree CherenkovTree("CherenkovTree", "");
    double CherenkovAngle_Reco_TrueEmissionPoint[2000], CherenkovAngle_Reco[2000],
           CherenkovAngle_True[2000], PhotonEnergy[2000];
    double OuterTracker_x, OuterTracker_y, OuterTracker_z;
    double BeforeRadiator_x, BeforeRadiator_y, BeforeRadiator_z;
    double Entrance_x, Entrance_y, Entrance_z;
    double MirrorHit_x, MirrorHit_y, MirrorHit_z;
    double RadiatorPosition_x, RadiatorPosition_y, RadiatorPosition_z;
    double CosTheta, Phi;
    int NumberPhotons = 0;
    std::size_t RadiatorRowNumber, RadiatorColumnNumber, TrackNumber;
    std::size_t FinalRadiatorRowNumber, FinalRadiatorColumnNumber;
    ParticleTrack::Location ParticleLocation;
    Photon::Status PhotonStatus[2000];
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
    CherenkovTree.Branch("PhotonEnergy",
			 &PhotonEnergy,
			 "PhotonEnergy[NumberPhotons]/D");
    CherenkovTree.Branch("FinalRadiatorColumnNumber",
			 &FinalRadiatorColumnNumber,
			 "FinalRadiatorColumnNumber/l");
    CherenkovTree.Branch("FinalRadiatorRowNumber",
			 &FinalRadiatorRowNumber,
			 "FinalRadiatorRowNumber/l");
    CherenkovTree.Branch("RadiatorRowNumber",
			 &RadiatorRowNumber,
			 "RadiatorRowNumber/l");
    CherenkovTree.Branch("RadiatorColumnNumber",
			 &RadiatorColumnNumber,
			 "RadiatorColumnNumber/l");
    CherenkovTree.Branch("TrackNumber",
			 &TrackNumber,
			 "TrackNumber/l");
    CherenkovTree.Branch("ParticleLocation", &ParticleLocation, "ParticleLocation/I");
    CherenkovTree.Branch("PhotonStatus", &PhotonStatus, "PhotonStatus[NumberPhotons]/I");
    CherenkovTree.Branch("Entrance_x", &Entrance_x);
    CherenkovTree.Branch("Entrance_y", &Entrance_y);
    CherenkovTree.Branch("Entrance_z", &Entrance_z);
    CherenkovTree.Branch("MirrorHit_x", &MirrorHit_x);
    CherenkovTree.Branch("MirrorHit_y", &MirrorHit_y);
    CherenkovTree.Branch("MirrorHit_z", &MirrorHit_z);
    CherenkovTree.Branch("OuterTracker_x", &OuterTracker_x);
    CherenkovTree.Branch("OuterTracker_y", &OuterTracker_y);
    CherenkovTree.Branch("OuterTracker_z", &OuterTracker_z);
    CherenkovTree.Branch("BeforeRadiator_x", &BeforeRadiator_x);
    CherenkovTree.Branch("BeforeRadiator_y", &BeforeRadiator_y);
    CherenkovTree.Branch("BeforeRadiator_z", &BeforeRadiator_z);
    CherenkovTree.Branch("RadiatorPosition_x", &RadiatorPosition_x);
    CherenkovTree.Branch("RadiatorPosition_y", &RadiatorPosition_y);
    CherenkovTree.Branch("RadiatorPosition_z", &RadiatorPosition_z);
    CherenkovTree.Branch("CosTheta", &CosTheta);
    CherenkovTree.Branch("Phi", &Phi);
    const std::size_t NumberTracks = Settings::GetSizeT("General/NumberTracks");
    const bool DrawAllTracks = Settings::GetBool("General/DrawAllTracks");
    const std::vector<int> TracksToDraw = Settings::GetIntVector("General/TrackToDraw");
    const std::size_t RowToDraw = Settings::GetSizeT("EventDisplay/RowToDraw");
    const bool DrawMissPhoton = Settings::GetBool("General/DrawMissPhoton");
    const bool Aerogel = Settings::GetString("General/GasOrAerogel") == "Aerogel";
    for(std::size_t i = 0; i < NumberTracks; i++) {
      NumberPhotons = 0;
      TrackNumber = i;
      RadiatorRowNumber = static_cast<std::size_t>(-1);
      RadiatorColumnNumber = static_cast<std::size_t>(-1);
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
      auto OuterTrackerPosition = particleTrack.GetPosition();
      OuterTracker_x = OuterTrackerPosition.X();
      OuterTracker_y = OuterTrackerPosition.Y();
      OuterTracker_z = OuterTrackerPosition.Z();
      if(!particleTrack.FindRadiator(*radiatorArray)) {
	continue;
      }
      auto BeforeRadiatorPosition = particleTrack.GetPosition();
      BeforeRadiator_x = BeforeRadiatorPosition.X();
      BeforeRadiator_y = BeforeRadiatorPosition.Y();
      BeforeRadiator_z = BeforeRadiatorPosition.Z();
      Phi = particleTrack.GetPosition().Phi();
      particleTrack.ConvertToRadiatorCoordinates();
      RadiatorRowNumber = particleTrack.GetRadiatorRowNumber();
      RadiatorColumnNumber = particleTrack.GetRadiatorColumnNumber();
      auto EntranceWindowPosition = particleTrack.GetEntranceWindowPosition();
      Entrance_x = EntranceWindowPosition.X();
      Entrance_y = EntranceWindowPosition.Y();
      Entrance_z = EntranceWindowPosition.Z();
      if(particleTrack.GetParticleLocation() != ParticleTrack::Location::EntranceWindow) {
	ParticleLocation = particleTrack.GetParticleLocation();
	CherenkovTree.Fill();
	continue;
      }
      particleTrack.TrackThroughAerogel();
      std::vector<Photon> Photons;
      if(Aerogel) {
	Photons = particleTrack.GeneratePhotonsFromAerogel();
      }
      particleTrack.TrackThroughGasToMirror();
      ParticleLocation = particleTrack.GetParticleLocation();
      if(particleTrack.GetParticleLocation() != ParticleTrack::Location::Mirror) {
	const bool HitMirror = particleTrack.TrackToNextCell(*radiatorArray);
	if(!HitMirror) {
	  continue;
	}
      }
      Phi = particleTrack.GetPosition().Phi();
      auto MirrorHitPosition = particleTrack.GetPosition();
      MirrorHit_x = MirrorHitPosition.X();
      MirrorHit_y = MirrorHitPosition.Y();
      MirrorHit_z = MirrorHitPosition.Z();
      auto RadiatorPosition = particleTrack.GetRadiatorCell()->GetRadiatorPosition();
      RadiatorPosition_x = RadiatorPosition.X();
      RadiatorPosition_y = RadiatorPosition.Y();
      RadiatorPosition_z = RadiatorPosition.Z();
      FinalRadiatorRowNumber = particleTrack.GetRadiatorRowNumber();
      FinalRadiatorColumnNumber = particleTrack.GetRadiatorColumnNumber();
      auto IsTrackDraw = [&] () {
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
      if(DrawThisTrack && FinalRadiatorRowNumber == RowToDraw) {
	eventDisplay.AddObject(particleTrack.DrawParticleTrack());
      }
      if(!Aerogel) {
	Photons = particleTrack.GeneratePhotonsFromGas();
      }
      for(auto &Photon : Photons) {
	if(NumberPhotons >= 2000) {
	  std::cout << "Warning! Number of photons is greater than 2000\n";
	}
	CherenkovAngle_True[NumberPhotons] = TMath::ACos(Photon.GetCosCherenkovAngle());
	auto photonHit = PhotonMapper::TracePhoton(Photon, *radiatorArray);
	if(DrawThisTrack && Photon.GetRadiatorRowNumber() == RowToDraw) {
	  eventDisplay.AddObject(Photon.DrawPhotonPath());
	}
	if(!Photon.GetMirrorHitPosition()) {
	  if(DrawMissPhoton && Photon.GetStatus() == Photon::Status::MirrorMiss) {
	    eventDisplay.AddObject(Photon.DrawPhotonPath());
	  }
	  CherenkovAngle_Reco_TrueEmissionPoint[NumberPhotons] = -1.0;
	  CherenkovAngle_Reco[NumberPhotons] = -1.0;
	} else {
	  auto Radiator = Aerogel ?
	                  Photon::Radiator::Aerogel :
	                  Photon::Radiator::Gas;
	  auto reconstructedPhoton =
	    PhotonReconstructor::ReconstructPhoton(particleTrack,
						   *photonHit,
						   Radiator);
	  CherenkovAngle_Reco_TrueEmissionPoint[NumberPhotons] =
	    TMath::ACos(reconstructedPhoton.m_CosCherenkovAngle_TrueEmissionPoint);
	  CherenkovAngle_Reco[NumberPhotons] =
	    TMath::ACos(reconstructedPhoton.m_CosCherenkovAngle);
	}
	PhotonEnergy[NumberPhotons] = Photon.GetEnergy();
	PhotonStatus[NumberPhotons] = Photon.GetStatus();
	NumberPhotons++;
      }
      CherenkovTree.Fill();
    }
    std::string EventDisplayFilename("EventDisplay");
    if(Settings::GetString("General/BarrelOrEndcap") != "Barrel") {
      EventDisplayFilename += "_EndCap";
    } else {
      if(Settings::GetInt("EventDisplay/RowToDraw") == 1) {
	EventDisplayFilename += "_MainRow";
      } else if(Settings::GetInt("EventDisplay/RowToDraw") == 2) {
	EventDisplayFilename += "_UpperRow";
      }
    }
    EventDisplayFilename += ".pdf";
    eventDisplay.DrawEventDisplay(EventDisplayFilename);
    CherenkovTree.Write();
    CherenkovFile.Close();
  }
  return 0;
}
