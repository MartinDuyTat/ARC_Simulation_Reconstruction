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
#include"ARCVector.h"

using Vector = ROOT::Math::XYZVector;

struct CherenkovFile {
  /**
   * The file
   */
  TFile File;
  /**
   * The TTree with all information
   */
  TTree CherenkovTree;
  /**
   * All the variables that are saved
   */
  double CherenkovAngle_Reco_TrueEmissionPoint[2000], CherenkovAngle_Reco[2000],
         CherenkovAngle_True[2000], PhotonEnergy[2000];
  double Entrance_x, Entrance_y, Entrance_z;
  double MirrorHit_x, MirrorHit_y, MirrorHit_z;
  double PhotonHit_x[2000], PhotonHit_y[2000], PhotonHit_z[2000];
  double Momentum, CosTheta, Phi;
  int NumberPhotons = 0;
  double NumberGoodPhotons = 0.0;
  int HasMigrated[2000];
  std::size_t RadiatorRowNumber, RadiatorColumnNumber, TrackNumber;
  std::size_t FinalRadiatorRowNumber, FinalRadiatorColumnNumber;
  ParticleTrack::Location ParticleLocation;
  Photon::Status PhotonStatus[2000];
  double SinglePhotonResolution, TotalResolution, Significance;
  /**
   * Constructor that prepares the file and tree
   */
  CherenkovFile(const std::string &Filename):
    File(Filename.c_str(), "RECREATE"),
    CherenkovTree("CherenkovTree", "") {
    CherenkovTree.Branch("NumberPhotons", &NumberPhotons);
    CherenkovTree.Branch("NumberGoodPhotons", &NumberGoodPhotons);
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
    CherenkovTree.Branch("HasMigrated",
			 &HasMigrated,
			 "HasMigrated[NumberPhotons]/I");
    CherenkovTree.Branch("ParticleLocation", &ParticleLocation, "ParticleLocation/I");
    CherenkovTree.Branch("PhotonStatus", &PhotonStatus, "PhotonStatus[NumberPhotons]/I");
    CherenkovTree.Branch("Entrance_x", &Entrance_x);
    CherenkovTree.Branch("Entrance_y", &Entrance_y);
    CherenkovTree.Branch("Entrance_z", &Entrance_z);
    CherenkovTree.Branch("MirrorHit_x", &MirrorHit_x);
    CherenkovTree.Branch("MirrorHit_y", &MirrorHit_y);
    CherenkovTree.Branch("MirrorHit_z", &MirrorHit_z);
    CherenkovTree.Branch("PhotonHit_x",
			 &PhotonHit_x,
			 "PhotonHit_x[NumberPhotons]/D");
    CherenkovTree.Branch("PhotonHit_y",
			 &PhotonHit_y,
			 "PhotonHit_y[NumberPhotons]/D");
    CherenkovTree.Branch("PhotonHit_z",
			 &PhotonHit_z,
			 "PhotonHit_z[NumberPhotons]/D");
    CherenkovTree.Branch("Momentum", &Momentum);
    CherenkovTree.Branch("CosTheta", &CosTheta);
    CherenkovTree.Branch("Phi", &Phi);
    CherenkovTree.Branch("SinglePhotonResolution", &SinglePhotonResolution);
    CherenkovTree.Branch("TotalResolution", &TotalResolution);
    CherenkovTree.Branch("Significance", &Significance);
  }
  /**
   * Fill the TTree
   */
  void Fill() {
    CherenkovTree.Fill();
  }
  /**
   * Save and close
   */
  void Close() {
    CherenkovTree.Write();
    File.Close();
  }
};

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
    ParticleTrack particleTrack(ParticleID, Momentum,
				0, InnerTracker.GetFieldStrength());
    if(!particleTrack.TrackThroughTracker(InnerTracker)) {
      std::cout << "Could not track through tracker...\n";
      return 0;
    }
    particleTrack.FindRadiator(*radiatorArray);
    particleTrack.ConvertToRadiatorCoordinates();
    if(!particleTrack.TrackThroughAerogel()) {
      std::cout << "Could not track through aerogel...\n";
      return 0;
    }
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
    CherenkovFile File("CherenkovFile.root");
    const std::size_t NumberTracks = Settings::GetSizeT("General/NumberTracks");
    const bool Aerogel = Settings::GetString("General/GasOrAerogel") == "Aerogel";
    const Photon::Radiator Radiator = Aerogel ?
                                      Photon::Radiator::Aerogel :
                                      Photon::Radiator::Gas;
    const int Hypothesis1 = Settings::GetInt("General/MassHypothesis1");
    const int Hypothesis2 = Settings::GetInt("General/MassHypothesis2");
    for(std::size_t i = 0; i < NumberTracks; i++) {
      File.NumberPhotons = 0;
      File.NumberGoodPhotons = 0.0;
      File.TrackNumber = i;
      File.RadiatorRowNumber = static_cast<std::size_t>(-1);
      File.RadiatorColumnNumber = static_cast<std::size_t>(-1);
      auto GetMomentum = [&] () {
	if(BarrelOrEndcap == "Barrel") {
	  return Utilities::GenerateRandomBarrelTrack(File.CosTheta, File.Phi);
	} else if(BarrelOrEndcap == "EndCap") {
	  return Utilities::GenerateRandomEndCapTrack();
	} else {
	  return Vector(0.0, 0.0, 0.0);
	}
      };
      const Vector Momentum = GetMomentum();
      File.Momentum = TMath::Sqrt(Momentum.Mag2());
      if(BarrelOrEndcap == "EndCap") {
	File.CosTheta = TMath::Cos(Momentum.Theta());
      }
      const Vector Position(0.0, 0.0, 0.0);
      const int ParticleID = Settings::GetInt("Particle/ID");
      ParticleTrack particleTrack(ParticleID, Momentum,
				  i, InnerTracker.GetFieldStrength());
      const double MomentumMag = TMath::Sqrt(Momentum.Mag2());
      const double CherenkovAngleDifference =
	Utilities::GetCherenkovAngleDifference(MomentumMag,
					       Hypothesis1,
					       Hypothesis2,
					       Radiator);
      if(!particleTrack.TrackThroughTracker(InnerTracker)) {
	continue;
      }
      if(!particleTrack.FindRadiator(*radiatorArray)) {
	continue;
      }
      File.Phi = particleTrack.GetPosition().GlobalVector().Phi();
      File.RadiatorRowNumber = particleTrack.GetRadiatorRowNumber();
      File.RadiatorColumnNumber = particleTrack.GetRadiatorColumnNumber();
      auto EntranceWindowPosition = particleTrack.GetEntranceWindowPosition();
      File.Entrance_x = EntranceWindowPosition.X();
      File.Entrance_y = EntranceWindowPosition.Y();
      File.Entrance_z = EntranceWindowPosition.Z();
      if(particleTrack.GetParticleLocation() != ParticleTrack::Location::EntranceWindow) {
	File.ParticleLocation = particleTrack.GetParticleLocation();
	File.Fill();
	continue;
      }
      if(!particleTrack.TrackThroughAerogel()) {
	continue;
      }
      std::vector<Photon> Photons;
      if(Aerogel) {
	Photons = particleTrack.GeneratePhotonsFromAerogel();
      }
      particleTrack.TrackThroughGasToMirror();
      if(particleTrack.GetParticleLocation() != ParticleTrack::Location::Mirror) {
	const bool HitMirror = particleTrack.TrackToNextCell(*radiatorArray);
	if(!HitMirror) {
	  File.ParticleLocation = particleTrack.GetParticleLocation();
	  File.FinalRadiatorRowNumber = static_cast<std::size_t>(-1);
	  File.FinalRadiatorColumnNumber = static_cast<std::size_t>(-1);
	  File.Fill();
	  continue;
	}
      }
      if(!Aerogel) {
	Photons = particleTrack.GeneratePhotonsFromGas();
      }
      File.ParticleLocation = particleTrack.GetParticleLocation();
      File.Phi = particleTrack.GetPosition().GlobalVector().Phi();
      auto MirrorHitPosition = particleTrack.GetPosition().GlobalVector();
      File.MirrorHit_x = MirrorHitPosition.X();
      File.MirrorHit_y = MirrorHitPosition.Y();
      File.MirrorHit_z = MirrorHitPosition.Z();
      File.FinalRadiatorRowNumber = particleTrack.GetRadiatorRowNumber();
      File.FinalRadiatorColumnNumber = particleTrack.GetRadiatorColumnNumber();
      eventDisplay.AddObject(particleTrack.DrawParticleTrack());
      std::vector<double> GoodAngles;
      for(auto &Photon : Photons) {
	if(File.NumberPhotons >= 2000) {
	  std::cout << "Warning! Number of photons is greater than 2000\n";
	}
	File.CherenkovAngle_True[File.NumberPhotons] =
	  TMath::ACos(Photon.GetCosCherenkovAngle());
	auto photonHit = PhotonMapper::TracePhoton(Photon, *radiatorArray);
	eventDisplay.AddObject(Photon.DrawPhotonPath());
	if(!Photon.GetMirrorHitPosition()) {
	  File.CherenkovAngle_Reco_TrueEmissionPoint[File.NumberPhotons] = -1.0;
	  File.CherenkovAngle_Reco[File.NumberPhotons] = -1.0;
	} else {
	  auto reconstructedPhoton =
	    PhotonReconstructor::ReconstructPhoton(*photonHit);
	  File.CherenkovAngle_Reco_TrueEmissionPoint[File.NumberPhotons] =
	    TMath::ACos(reconstructedPhoton.m_CosCherenkovAngle_TrueEmissionPoint);
	  File.CherenkovAngle_Reco[File.NumberPhotons] =
	    TMath::ACos(reconstructedPhoton.m_CosCherenkovAngle);
	  File.HasMigrated[File.NumberPhotons] = Photon.HasPhotonMigrated() ? 1 : 0;
	  if(Photon.GetStatus() == Photon::Status::DetectorHit) {
	    File.NumberGoodPhotons += Photon.GetWeight();
	    GoodAngles.push_back(File.CherenkovAngle_Reco[File.NumberPhotons]);
	  }
	}
	File.PhotonEnergy[File.NumberPhotons] = Photon.GetEnergy();
	File.PhotonStatus[File.NumberPhotons] = Photon.GetStatus();
	const auto PhotonPosition = Photon.GetPosition().GlobalVector();
	File.PhotonHit_x[File.NumberPhotons] = PhotonPosition.X();
	File.PhotonHit_y[File.NumberPhotons] = PhotonPosition.Y();
	File.PhotonHit_z[File.NumberPhotons] = PhotonPosition.Z();
	File.NumberPhotons++;
      }
      const double RMS = TMath::RMS(GoodAngles.begin(), GoodAngles.end());
      File.SinglePhotonResolution = RMS;
      File.TotalResolution = RMS/TMath::Sqrt(File.NumberGoodPhotons);
      File.Significance = CherenkovAngleDifference/File.TotalResolution;
      File.Fill();
    }
    eventDisplay.DrawEventDisplay();
    File.Close();
  }
  return 0;
}
