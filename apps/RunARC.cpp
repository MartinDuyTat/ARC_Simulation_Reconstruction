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

Vector VectorFromSpherical(double R, double CosTheta, double Phi);

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
						Settings::GetDouble("Particle/CosTheta"),
						Settings::GetDouble("Particle/Phi"));
    const int ParticleID = Settings::GetInt("Particle/ID");;
    ParticleTrack particleTrack(ParticleID, Momentum);
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
    radiatorArray(0, 0)->GetDetector().PlotHits("PhotonHits.pdf");
  } else if(RunMode == "CherenkovAngleResolution") {
    std::cout << "Run mode: Cherenkov angle resolution\n";
    TFile CherenkovFile("CherenkovFile.root", "RECREATE");
    TTree CherenkovTree("CherenkovTree", "");
    double CherenkovAngle_Reco_TrueEmissionPoint[200], CherenkovAngle_Reco[200],
           CherenkovAngle_True[200], PhotonEnergy[200];
    double Entrance_x, Entrance_y, Entrance_z;
    double CosTheta, Phi;
    int NumberPhotons = 0;
    int RadiatorRowNumber[200], RadiatorColumnNumber[200], TrackNumber;
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
    CherenkovTree.Branch("RadiatorRowNumber",
			 &RadiatorRowNumber,
			 "RadiatorRowNumber[NumberPhotons]/D");
    CherenkovTree.Branch("RadiatorColumnNumber",
			 &RadiatorColumnNumber,
			 "RadiatorColumnNumber[NumberPhotons]/D");
    CherenkovTree.Branch("TrackNumber", &TrackNumber);
    CherenkovTree.Branch("ParticleLocation", &ParticleLocation, "ParticleLocation/I");
    CherenkovTree.Branch("PhotonStatus", &PhotonStatus, "PhotonStatus[NumberPhotons]/I");
    CherenkovTree.Branch("Entrance_x", &Entrance_x);
    CherenkovTree.Branch("Entrance_y", &Entrance_y);
    CherenkovTree.Branch("Entrance_z", &Entrance_z);
    CherenkovTree.Branch("CosTheta", &CosTheta);
    CherenkovTree.Branch("Phi", &Phi);
    const int NumberTracks = Settings::GetInt("General/NumberTracks");
    const std::vector<int> TracksToDraw = Settings::GetIntVector("General/TrackToDraw");
    const bool DrawMissPhoton = Settings::GetBool("General/DrawMissPhoton");
    for(int i = 0; i < NumberTracks; i++) {
      NumberPhotons = 0;
      TrackNumber = i;
      const bool DrawThisTrack = std::find(TracksToDraw.begin(),
					   TracksToDraw.end(), i) != TracksToDraw.end();
      const double Radius = Settings::GetDouble("ARCGeometry/Radius");
      const double z = gRandom->Uniform(Settings::GetDouble("Particle/z_min"),
					Settings::GetDouble("Particle/z_max"));
      CosTheta = z/TMath::Sqrt(z*z + Radius*Radius);
      Phi = Settings::GetBool("Particle/RandomPhi")
	  ? gRandom->Uniform(-TMath::Pi(), TMath::Pi())
	  : gRandom->Uniform(Settings::GetDouble("Particle/Phi_min"),
	                     Settings::GetDouble("Particle/Phi_max"));
      const double MomentumMag = Settings::GetDouble("Particle/Momentum");
      const Vector Momentum = VectorFromSpherical(MomentumMag, CosTheta, Phi);
      const Vector Position(0.0, 0.0, 0.0);
      const int ParticleID = Settings::GetInt("Particle/ID");;
      ParticleTrack particleTrack(ParticleID, Momentum, Position);
      particleTrack.TrackThroughTracker(InnerTracker);
      particleTrack.ConvertToRadiatorCoordinates(radiatorArray);
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
	CherenkovAngle_True[NumberPhotons] = Photon.m_CherenkovAngle;
	PhotonMapper::TracePhoton(Photon);
	if(DrawThisTrack) {
	  eventDisplay.AddObject(Photon.DrawPhotonPath());
	}
	if(!Photon.m_MirrorHitPosition) {
	  if(DrawMissPhoton && Photon.m_Status == Photon::Status::MirrorMiss) {
	    eventDisplay.AddObject(Photon.DrawPhotonPath());
	  }
	  continue;
	}
	auto reconstructedPhoton =
	  PhotonReconstructor::ReconstructPhoton(particleTrack,
						 particleTrack.GetPhotonHits().back(),
						 Photon::Radiator::Gas);
	PhotonEnergy[NumberPhotons] = Photon.m_Energy;
	PhotonStatus[NumberPhotons] = Photon.m_Status;
	RadiatorRowNumber[NumberPhotons] = Photon.m_RadiatorCell->GetCellNumber().first;
	RadiatorColumnNumber[NumberPhotons] = Photon.m_RadiatorCell->GetCellNumber().second;
	CherenkovAngle_Reco_TrueEmissionPoint[NumberPhotons] =
	  reconstructedPhoton.m_CherenkovAngle_TrueEmissionPoint;
	CherenkovAngle_Reco[NumberPhotons] = reconstructedPhoton.m_CherenkovAngle;
	NumberPhotons++;
      }
      CherenkovTree.Fill();
    }
    eventDisplay.DrawEventDisplay("EventDisplay.pdf");
    CherenkovTree.Write();
    CherenkovFile.Close();
  }
  return 0;
}

Vector VectorFromSpherical(double R, double CosTheta, double Phi) {
  const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
  const double CosPhi = TMath::Cos(Phi);
  const double SinPhi = TMath::Sin(Phi);
  return Vector{R*CosPhi*SinTheta, R*SinPhi*SinTheta, R*CosTheta};
}
