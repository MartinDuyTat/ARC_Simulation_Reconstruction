// Martin Duy Tat 29th April 2022

#include<stdexcept>
#include"TMath.h"
#include"TRandom.h"
#include"Math/RotationY.h"
#include"Math/RotationZ.h"
#include"ParticleTrack.h"
#include"Photon.h"
#include"TrackingVolume.h"
#include"RadiatorCell.h"
#include"ParticleMass.h"
#include"Settings.h"
#include"SiPM.h"
#include"RadiatorArray.h"

ParticleTrack::ParticleTrack(const Vector &Momentum, int ParticleID): m_Momentum(Momentum),
								      m_Position(0.0, 0.0, 0.0),
								      m_InitialPosition(0.0, 0.0, 0.0),
								      m_ParticleID(ParticleID),
								      m_TrackedThroughTracker(false),
								      m_TrackedThroughRadiator(false),
								      m_AerogelEntry(0.0, 0.0, 0.0),
								      m_AerogelExit(0.0, 0.0, 0.0),
								      m_GasEntry(0.0, 0.0, 0.0),
								      m_GasExit(0.0, 0.0, 0.0),
								      m_CoordinateSystem(CoordinateSystem::GlobalDetector) {
}

void ParticleTrack::TrackThroughTracker(const TrackingVolume &InnerTracker) {
  // TODO: Account for magnetic field
  if(m_TrackedThroughTracker) {
    throw std::runtime_error("Cannot track particle through inner tracker again");
  }
  m_Position += m_Momentum.Unit()*InnerTracker.GetRadius()/TMath::Sin(m_Momentum.Theta());
  m_TrackedThroughTracker = true;
}

void ParticleTrack::ConvertToRadiatorCoordinates(RadiatorArray &radiatorArray) {
  // Check if coordinate system if correct
  if(m_CoordinateSystem == CoordinateSystem::LocalRadiator) {
    throw std::runtime_error("Particle position is already in local radiator coordinates");
  }
  // First rotate in azimuthal direction to map to cell near phi = 0
  MapPhiBack(m_Position);
  MapPhiBack(m_Momentum);
  // Then rotate around y-axis so that z axis now points towards the high pT cell
  RotateY(m_Position);
  RotateY(m_Momentum);
  // Figure out which radiator cell the track goes through
  m_RadiatorCell = radiatorArray.WhichRadiator(m_Position);
  // Finally shift coordinates so that the origin is the detector plane of the radiator cell
  m_Position -= m_RadiatorCell->GetRadiatorPosition();
  m_CoordinateSystem = CoordinateSystem::LocalRadiator;
  // Check if particle is within acceptance
  if(!m_RadiatorCell->IsInsideThetaBoundary(m_Position)) {
    throw std::runtime_error("Particle outside of radiator acceptance");
  }
}

void ParticleTrack::TrackThroughRadiatorCell() {
  // TODO: Allow particle to go through multiple cells
  if(m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot track particle through radiator cell again");
  }
  if(!m_TrackedThroughTracker) {
    throw std::runtime_error("Particle must be tracked through inner tracker before tracking through radiator cell");
  }
  const double ZDistToAerogel = -m_Position.Z();
  const double Slope = 1.0/TMath::Cos(m_Momentum.Theta());
  m_Position += m_Momentum.Unit()*Slope*ZDistToAerogel;
  m_AerogelEntry = m_Position;
  m_Position += m_Momentum.Unit()*Slope*m_RadiatorCell->GetAerogelThickness();
  m_AerogelExit = m_Position;
  m_GasEntry = m_AerogelExit;
  TrackThroughGasToMirror();
}

void ParticleTrack::TrackThroughGasToMirror() {
  do {
    // Solve quadratic s^2 - 2sb + c = 0 to find interserction of particle track and mirror
    auto MirrorCentre = m_RadiatorCell->GetMirrorCentre();
    const double MirrorRadius = m_RadiatorCell->GetMirrorCurvature();
    const auto Direction = m_Momentum.Unit();
    const double b = (MirrorCentre - m_Position).Dot(Direction);
    const double c = MirrorCentre.Mag2() + m_Position.Mag2() - MirrorRadius*MirrorRadius - 2*m_Position.Dot(MirrorCentre);
    const double s1 = b + TMath::Sqrt(b*b - c);
    const double s2 = b - TMath::Sqrt(b*b - c);
    // Pick forwards moving particle solution
    const double s = std::max(s1, s2);
    m_Position += Direction*s;
    m_GasExit = m_Position;
    m_TrackedThroughRadiator = m_RadiatorCell->IsInsideThetaBoundary(m_Position);
    if(!m_TrackedThroughRadiator) {
      if(!SwapRadiatorCell()) {
	m_TrackedThroughRadiator = true;
	break;
      }
    }
  } while(!m_TrackedThroughRadiator);
}

bool ParticleTrack::SwapRadiatorCell() {
  if(m_RadiatorCell->IsInsideThetaBoundary(m_Position)) {
    return true;
  }
  if(m_RadiatorCell->GetFirstMiddleLast() == RadiatorCell::FirstMiddleLast::Single) {
    return false;
  }
  auto CurrentRadiatorPosition = m_RadiatorCell->GetRadiatorPosition();
  if(m_Position.X() > m_RadiatorCell->GetThetaLength()/2.0) {
    if(m_RadiatorCell->GetFirstMiddleLast() == RadiatorCell::FirstMiddleLast::Last) {
      return false;
    }
    m_RadiatorCell++;
  } else if(m_Position.X() < -m_RadiatorCell->GetThetaLength()/2.0) {
    if(m_RadiatorCell->GetFirstMiddleLast() == RadiatorCell::FirstMiddleLast::First) {
      return false;
    }
    m_RadiatorCell--;
  } else {
    throw std::runtime_error("SwapRadiatorCell cannot figure out where particle track is");
  }
  auto NewRadiatorPosition = m_RadiatorCell->GetRadiatorPosition();
  ChangeCoordinateOrigin(NewRadiatorPosition - CurrentRadiatorPosition);
  return true;
}

void ParticleTrack::ChangeCoordinateOrigin(const Vector &Shift) {
  m_Position -= Shift;
  m_AerogelEntry -= Shift;
  m_AerogelExit -= Shift;
  m_GasEntry -= Shift;
  m_GasExit -= Shift;
}

Photon ParticleTrack::GeneratePhotonFromAerogel() const {
  if(!m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot generate photons from tracks that have not been tracked through the radiator");
  }
  Photon photon = GeneratePhoton(m_AerogelEntry, m_AerogelExit, Photon::Radiator::Aerogel);
  return photon;
}

Photon ParticleTrack::GeneratePhotonFromGas() const {
  if(!m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot generate photons from tracks that have not been tracked through the radiator");
  }
  Photon photon = GeneratePhoton(m_GasEntry, m_GasExit, Photon::Radiator::Gas);
  return photon;
}

std::vector<Photon> ParticleTrack::GeneratePhotonsFromAerogel() const {
  if(!m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot generate photons from tracks that have not been tracked through the radiator");
  }
  const double RadiatorDistance = TMath::Sqrt((m_AerogelExit - m_AerogelEntry).Mag2());
  const int PhotonYield = std::round(GetPhotonYield(RadiatorDistance, Beta(), 1.03));
  std::vector<Photon> Photons;
  for(int i = 0; i < PhotonYield; i++) {
    Photons.push_back(GeneratePhoton(m_AerogelEntry, m_AerogelExit, Photon::Radiator::Aerogel));
  }
  return Photons;
}

std::vector<Photon> ParticleTrack::GeneratePhotonsFromGas() const {
  if(!m_TrackedThroughRadiator) {
    throw std::runtime_error("Cannot generate photons from tracks that have not been tracked through the radiator");
  }
  const double RadiatorDistance = TMath::Sqrt((m_GasExit - m_GasEntry).Mag2());
  const int PhotonYield = std::round(GetPhotonYield(RadiatorDistance, Beta(), 1.0049));
  std::vector<Photon> Photons;
  for(int i = 0; i < PhotonYield; i++) {
    Photons.push_back(GeneratePhoton(m_GasEntry, m_GasExit, Photon::Radiator::Gas));
  }
  return Photons;
}

double ParticleTrack::GetIndexRefraction(Photon::Radiator Radiator, double Energy) const {
  switch(Radiator) {
    case Photon::Radiator::Aerogel:
      return 1.03;
    case Photon::Radiator::Gas:
      {
	const double Lambda = 1239.841987427/Energy;
	if(Settings::GetBool("General/ChromaticDispersion")) {
	  // Pressure at 3.5 bar
	  // Sellmeier equation with coefficients from https://twiki.cern.ch/twiki/bin/view/LHCb/C4F10
	  // They are similar to A. Filippas, et al. Nucl. Instr. and Meth. B, 196 (2002), p. 340 but now quite...?
	  return 3.5*0.25324*1e-6/((1.0/(73.7*73.7)) - (1.0/(Lambda*Lambda)));
	} else {
	  return 1.0049;
	}
      }
    default:
      return 1.0;
  }
}

Photon ParticleTrack::GeneratePhoton(const Vector &Entry, const Vector &Exit, Photon::Radiator Radiator) const {
  // TODO: Account for dispersion
  const double Energy = gRandom->Uniform(1.78, 4.31);
  const double n_phase = GetIndexRefraction(Radiator, Energy);
  const double RandomFraction = Settings::GetBool("General/RandomEmissionPoint") ? gRandom->Uniform(0.005, 0.995) : 0.5;
  const Vector EmissionPoint = Entry + (Exit - Entry)*RandomFraction;
  const double phi = gRandom->Uniform(0.0, 2*TMath::Pi());
  const double CosTheta = 1.0/(Beta()*n_phase);
  const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
  Vector Direction(SinTheta*TMath::Cos(phi), SinTheta*TMath::Sin(phi), CosTheta);
  // Rotate to particle frame
  const ROOT::Math::RotationY RotateY(m_Momentum.Theta());
  Direction = RotateY(Direction);
  const ROOT::Math::RotationZ RotateZ(m_Momentum.Phi());
  Direction = RotateZ(Direction);
  return {EmissionPoint, Direction, Energy, TMath::ACos(CosTheta), Radiator, &(*m_RadiatorCell)};
}

double ParticleTrack::Beta() const {
  const double Momentum = TMath::Sqrt(m_Momentum.Mag2());
  const double Mass = ParticleMass::GetMass(m_ParticleID);
  return Momentum/TMath::Sqrt(Mass*Mass + Momentum*Momentum);
}

const Vector& ParticleTrack::GetMomentum() const {
  return m_Momentum;
}

const Vector& ParticleTrack::GetEntryPoint(Photon::Radiator Radiator) const {
  if(Radiator == Photon::Radiator::Gas) {
    return m_GasEntry;
  } else if(Radiator == Photon::Radiator::Aerogel) {
    return m_AerogelEntry;
  } else {
    throw std::runtime_error("Cannot find entry point of unknown radiator");
  }
}

const Vector& ParticleTrack::GetExitPoint(Photon::Radiator Radiator) const {
  if(Radiator == Photon::Radiator::Gas) {
    return m_GasExit;
  } else if(Radiator == Photon::Radiator::Aerogel) {
    return m_AerogelExit;
  } else {
    throw std::runtime_error("Cannot find entry point of unknown radiator");
  }
}

double ParticleTrack::GetPhotonYield(double x, double Beta, double n) const {
  // TODO: Move this to separate class
  const double Efficiency = 0.432;
  const double DeltaE = 2.55;
  return x*DeltaE*37000.0*(1.0 - 1.0/TMath::Power(Beta*n, 2))*Efficiency;
}

void ParticleTrack::MapPhiBack(Vector &Vec) const {
  const int PhiCells = Settings::GetInt("ARCGeometry/PhiCells");
  const double DeltaPhi = 2.0*TMath::Pi()/PhiCells;
  const int Sign = Vec.Phi() > 0 ? 1 : -1;
  const int PhiUnits = Sign*static_cast<int>((Sign*Vec.Phi() + 0.5*DeltaPhi)/DeltaPhi);
  const ROOT::Math::RotationZ RotateZ(-PhiUnits*DeltaPhi);
  Vec = RotateZ(Vec);
}

void ParticleTrack::RotateY(Vector &Vec) const {
  const double Temp = Vec.X();
  Vec.SetX(-Vec.Z());
  Vec.SetZ(Temp);
}

std::unique_ptr<TLine> ParticleTrack::DrawParticleTrack() const {
  const auto CurrentPosition = m_RadiatorCell->GetRadiatorPosition() + m_Position;
  TLine Track(m_InitialPosition.X(), m_InitialPosition.Z(), CurrentPosition.X(), CurrentPosition.Z());
  Track.SetLineColor(kRed);
  return std::make_unique<TLine>(Track);
}

const std::vector<PhotonHit>& ParticleTrack::GetPhotonHits() const {
  return m_RadiatorCell->GetDetector().GetPhotonHits();
}

const Vector& ParticleTrack::GetPosition() const {
  return m_Position;
}
