// Martin Duy Tat 29th April 2022

#include<stdexcept>
#include"TMath.h"
#include"TRandom.h"
#include"ParticleTrack.h"
#include"Photon.h"
#include"TrackingVolume.h"
#include"RadiatorCell.h"
#include"ParticleMass.h"

ParticleTrack::ParticleTrack(const Vector &Momentum, int ParticleID): m_Momentum(Momentum),
								      m_Position(0.0, 0.0, 0.0),
								      m_ParticleID(ParticleID),
								      m_TrackedThroughTracker(false),
								      m_TrackedThroughRadiator(false),
								      m_AerogelEntry(0.0, 0.0, 0.0),
								      m_AerogelExit(0.0, 0.0, 0.0),
								      m_GasEntry(0.0, 0.0, 0.0),
								      m_GasExit(0.0, 0.0, 0.0),
								      m_CoordinateSystem(CoordinateSystem::GlobalDetector) {
  gRandom->SetSeed(42);
}

void ParticleTrack::TrackThroughTracker(const TrackingVolume &InnerTracker) {
  // TODO: Account for magnetic field
  if(m_TrackedThroughTracker) {
    throw std::runtime_error("Cannot track particle through inner tracker again");
  }
  m_Position += m_Momentum.Unit()*InnerTracker.GetRadius();
  m_TrackedThroughTracker = true;
}

void ParticleTrack::ConvertToRadiatorCoordinates(const RadiatorCell &Cell) {
  if(m_CoordinateSystem == CoordinateSystem::LocalRadiator) {
    throw std::runtime_error("Particle position is already in local radiator coordinates");
  }
  m_Position -= Cell.GetRadiatorPosition();
  m_CoordinateSystem = CoordinateSystem::LocalRadiator;
}

void ParticleTrack::TrackThroughRadiatorCell(const RadiatorCell &Cell) {
  // TODO: Allow for tracks in all directions
  if(m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot track particle through radiator cell again");
  }
  if(!m_TrackedThroughTracker) {
    throw std::runtime_error("Particle must be tracked through inner tracker before tracking through radiator cell");
  }
  m_AerogelEntry = Vector(0.0, 0.0, 0.0);
  m_AerogelExit = Vector(0.0, 0.0, Cell.GetAerogelThickness());
  m_GasEntry = m_AerogelExit;
  double GasThickness = Cell.GetRadiatorThickness() - 2*Cell.GetVesselThickness() - Cell.GetCoolingThickness() - Cell.GetAerogelThickness();
  m_GasExit = m_GasEntry + Vector(0.0, 0.0, GasThickness);
  m_TrackedThroughRadiator = true;
}

Photon ParticleTrack::GeneratePhotonFromAerogel() const {
  if(!m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot generate photons from tracks that have not been tracked through the radiator");
  }
  Photon photon = GeneratePhoton(m_AerogelEntry, m_AerogelExit, 1.03);
  photon.m_Radiator = Photon::Radiator::Aerogel;
  return photon;
}

Photon ParticleTrack::GeneratePhotonFromGas() const {
  if(!m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot generate photons from tracks that have not been tracked through the radiator");
  }
  Photon photon = GeneratePhoton(m_GasEntry, m_GasExit, 1.0049);
  photon.m_Radiator = Photon::Radiator::Gas;
  return photon;
}

std::vector<Photon> ParticleTrack::GeneratePhotonsFromAerogel() const {
  if(!m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot generate photons from tracks that have not been tracked through the radiator");
  }
  const double RadiatorDistance = TMath::Sqrt((m_AerogelExit - m_AerogelEntry).Mag2());
  const int PhotonYield = GetPhotonYield(RadiatorDistance, Beta(), 1.03);
  std::vector<Photon> Photons;
  for(int i = 0; i < PhotonYield; i++) {
    Photons.push_back(GeneratePhoton(m_AerogelEntry, m_AerogelExit, 1.05));
    Photons.back().m_Radiator = Photon::Radiator::Aerogel;
  }
  return Photons;
}

std::vector<Photon> ParticleTrack::GeneratePhotonsFromGas() const {
  if(!m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot generate photons from tracks that have not been tracked through the radiator");
  }
  const double RadiatorDistance = TMath::Sqrt((m_GasExit - m_GasEntry).Mag2());
  const int PhotonYield = GetPhotonYield(RadiatorDistance, Beta(), 1.0049);
  std::vector<Photon> Photons;
  for(int i = 0; i < PhotonYield; i++) {
    Photons.push_back(GeneratePhoton(m_GasEntry, m_GasExit, 1.0049));
    Photons.back().m_Radiator = Photon::Radiator::Gas;
  }
  return Photons;
}

Photon ParticleTrack::GeneratePhoton(const Vector &Entry, const Vector &Exit, double n_phase) const {
  // TODO: Rotate between local particle coordinate, beamline coordinates and local radiator coordinates
  // TODO: Account for dispersion
  const double Length = TMath::Sqrt((Exit - Entry).Mag2());
  const double RandomFraction = gRandom->Uniform(0, 1);
  const Vector EmissionPoint = Entry + (Exit - Entry)*Length*RandomFraction;
  const double phi = gRandom->Uniform(0.0, 2*TMath::Pi());
  const double CosTheta = 1.0/(Beta()*n_phase);
  const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
  const Vector Direction(SinTheta*TMath::Cos(phi), SinTheta*TMath::Sin(phi), CosTheta);
  return Photon{EmissionPoint, Direction, 0.0, TMath::ACos(CosTheta)};
}

double ParticleTrack::Beta() const {
  const double Momentum = TMath::Sqrt(m_Momentum.Mag2());
  const double Mass = ParticleMass::GetMass(m_ParticleID);
  return Momentum/TMath::Sqrt(Mass*Mass + Momentum*Momentum);
}

const Vector& ParticleTrack::GetMomentum() const {
  return m_Momentum;
}

double ParticleTrack::GetPhotonYield(double x, double Beta, double n) const {
  // TODO: Move this to separate class
  const double Efficiency = 0.20;
  const double DeltaE = 4.0;
  return x*DeltaE*37000.0*(1.0 - 1.0/TMath::Power(Beta*n, 2))*Efficiency;
}
