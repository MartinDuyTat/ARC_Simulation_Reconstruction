// Martin Duy Tat 1st May 2022
/**
 * RunARC is an application for running ARC simulations and reconstructions
 * Run with RunARC <run mode>
 * Possible run modes:
 * "SingleTrack" Send a single charged track and generate a photon hit map
 * "CherenkovAngleResolution" Generate photons from a single track, reconstruct the photons and study the Cherenkov angle resolution
 */

#include<string>
#include<vector>
#include"TFile.h"
#include"TTree.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"ParticleTrack.h"
#include"TrackingVolume.h"
#include"Photon.h"
#include"PhotonMapper.h"
#include"PhotonReconstructor.h"
#include"RadiatorCell.h"

using Vector = ROOT::Math::XYZVector;

int main(int argc, char *argv[]) {
  if(argc != 2) {
    return 0;
  }
  const std::string RunMode(argv[1]);
  const Vector Momentum(0.0, 0.0, 5.0);
  const int ParticleID = 211;
  ParticleTrack particleTrack(Momentum, ParticleID);
  const TrackingVolume InnerTracker(1.05, 2.0);
  RadiatorCell radiatorCell(Vector(0.0, 0.0, 1.09));
  particleTrack.TrackThroughTracker(InnerTracker);
  particleTrack.ConvertToRadiatorCoordinates(radiatorCell);
  particleTrack.TrackThroughRadiatorCell(radiatorCell);
  if(RunMode == "SingleTrack") {
    auto PhotonsAerogel = particleTrack.GeneratePhotonsFromAerogel();
    auto PhotonsGas = particleTrack.GeneratePhotonsFromGas();
    for(auto &photon : PhotonsAerogel) {
      PhotonMapper::TracePhoton(photon, radiatorCell);
    }
    for(auto &photon : PhotonsGas) {
      PhotonMapper::TracePhoton(photon, radiatorCell);
    }
    radiatorCell.m_Detector.PlotHits("PhotonHits.png");
  } else if(RunMode == "CherenkovAngleResolution") {
    TFile CherenkovFile("CherenkovFile.root", "RECREATE");
    TTree CherenkovTree("CherenkovTree", "");
    double CherenkovAngle_Reco, CherenkovAngle_True;
    CherenkovTree.Branch("CherenkovAngle_Reco", &CherenkovAngle_Reco);
    CherenkovTree.Branch("CherenkovAngle_True", &CherenkovAngle_True);
    std::vector<Photon> Photons;
    for(int i = 0; i < 1000; i++) {
      Photons.push_back(particleTrack.GeneratePhotonFromGas());
      CherenkovAngle_True = Photons.back().m_CherenkovAngle;
      PhotonMapper::TracePhoton(Photons.back(), radiatorCell);
      auto reconstructedPhoton = PhotonReconstructor::ReconstructPhoton(particleTrack, radiatorCell.m_Detector.GetPhotonHits().back(), radiatorCell);
      CherenkovAngle_Reco = reconstructedPhoton.m_CherenkovAngle_TrueEmission_TrueIndexRefraction;
      CherenkovTree.Fill();
    }
    CherenkovTree.Write();
    CherenkovFile.Close();
  }
  return 0;
}
