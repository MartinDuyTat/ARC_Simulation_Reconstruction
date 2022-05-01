// Martin Duy Tat 1st May 2022
/**
 * RunARC is an application for running ARC simulations and reconstructions
 */

#include<vector>
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"ParticleTrack.h"
#include"TrackingVolume.h"
#include"Photon.h"
#include"PhotonMapper.h"
#include"RadiatorCell.h"

using Vector = ROOT::Math::XYZVector;

int main() {
  const Vector Momentum(0.0, 0.0, 5.0);
  const int ParticleID = 211;
  ParticleTrack particleTrack(Momentum, ParticleID);
  const TrackingVolume InnerTracker(1.05, 2.0);
  RadiatorCell radiatorCell(Vector(0.0, 0.0, 1.09));
  particleTrack.TrackThroughTracker(InnerTracker);
  particleTrack.ConvertToRadiatorCoordinates(radiatorCell);
  particleTrack.TrackThroughRadiatorCell(radiatorCell);
  auto PhotonsAerogel = particleTrack.GeneratePhotonsFromAerogel();
  auto PhotonsGas = particleTrack.GeneratePhotonsFromGas();
  for(auto &photon : PhotonsAerogel) {
    PhotonMapper::TracePhoton(photon, radiatorCell);
  }
  for(auto &photon : PhotonsGas) {
    PhotonMapper::TracePhoton(photon, radiatorCell);
  }
  radiatorCell.m_Detector.PlotHits("PhotonHits.png");
  return 0;
}
