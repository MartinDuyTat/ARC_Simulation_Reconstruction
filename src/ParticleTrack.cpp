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

ParticleTrack::ParticleTrack(int ParticleID,
			     const Vector &Momentum,
			     const Vector &Position):
  m_Momentum(Momentum),
  m_Position(Position),
  m_InitialPosition(Position),
  m_ParticleID(ParticleID),
  m_Location(Location::TrackerVolume),
  m_CoordinateSystem(CoordinateSystem::GlobalDetector),
  m_RandomEmissionPoint(Settings::GetBool("General/RandomEmissionPoint")),
  m_ChromaticDispersion(Settings::GetBool("General/ChromaticDispersion")),
  m_Mass(ParticleMass::GetMass(m_ParticleID)) {
}

void ParticleTrack::TrackThroughTracker(const TrackingVolume &InnerTracker) {
  // TODO: Account for magnetic field
  if(m_Location != Location::TrackerVolume) {
    throw std::runtime_error("Cannot track particle through inner tracker again");
  }
  // Solve quadratic to track particle to edge of tracking barrel
  const double b = m_Position.X()*m_Momentum.Unit().X()
                 + m_Position.Y()*m_Momentum.Unit().Y();
  const double c = m_Position.X()*m_Position.X()
                 + m_Position.Y()*m_Position.Y()
                 - InnerTracker.GetRadius()*InnerTracker.GetRadius();
  const double s = TMath::Sqrt(b*b - c) - b;
  const double xyUnit = TMath::Sqrt(TMath::Power(m_Momentum.Unit().X(), 2) +
				    TMath::Power(m_Momentum.Unit().Y(), 2));
  m_Position += s*m_Momentum.Unit()/xyUnit;
  m_Location = Location::EntranceWindow;
}

void ParticleTrack::FindRadiator(RadiatorArray &radiatorArray) {
  m_RadiatorCell = &(*radiatorArray.FindRadiator(*this));
}

void ParticleTrack::SetRadiator(const RadiatorCell &radiatorCell) {
  m_RadiatorCell = &radiatorCell;
}

void ParticleTrack::ConvertToRadiatorCoordinates() {
  // Check if coordinate system if correct
  if(m_CoordinateSystem == CoordinateSystem::LocalRadiator) {
    throw std::runtime_error("Particle position is already in local radiator coordinates");
  }
  // Then rotate around y-axis so that z axis now points towards the high pT cell
  SwapXZ(m_InitialPosition);
  SwapXZ(m_Position);
  SwapXZ(m_Momentum);
  // Finally shift coordinates so that the origin is the detector plane of the radiator cell
  m_Position -= m_RadiatorCell->GetRadiatorPosition();
  m_CoordinateSystem = CoordinateSystem::LocalRadiator;
  // Check if particle is within acceptance
  if(!m_RadiatorCell->IsInsideCell(m_Position)) {
    m_Location = Location::MissedEntranceWindow;
  }
  // Save the entrance window position
  m_EntranceWindowPosition = m_Position;
}

void ParticleTrack::TrackThroughRadiatorCell() {
  if(m_Location != Location::EntranceWindow) {
    throw std::runtime_error("Particle not at entrance window");
  }
  const double ZDistToAerogel = -m_Position.Z();
  const double Slope = 1.0/TMath::Cos(m_Momentum.Theta());
  m_Position += m_Momentum.Unit()*Slope*ZDistToAerogel;
  m_AerogelEntry = m_Position;
  m_Position += m_Momentum.Unit()*Slope*m_RadiatorCell->GetAerogelThickness();
  m_AerogelExit = m_Position;
  m_GasEntry = m_AerogelExit;
  m_Location = m_RadiatorCell->IsInsideCell(m_Position) ?
               Location::Radiator : Location::MissedRadiator;
  TrackThroughGasToMirror();
}

void ParticleTrack::TrackThroughGasToMirror() {
  if(m_Location != Location::Radiator) {
    return;
  }
  // Solve quadratic s^2 - 2sb + c = 0 to find interserction of particle track and mirror
  auto MirrorCentre = m_RadiatorCell->GetMirrorCentre();
  const double MirrorRadius = m_RadiatorCell->GetMirrorCurvature();
  const auto Direction = m_Momentum.Unit();
  const double b = (MirrorCentre - m_Position).Dot(Direction);
  const double c = MirrorCentre.Mag2()
                 + m_Position.Mag2()
                 - MirrorRadius*MirrorRadius
                 - 2*m_Position.Dot(MirrorCentre);
  const double s1 = b + TMath::Sqrt(b*b - c);
  const double s2 = b - TMath::Sqrt(b*b - c);
  // Pick forwards moving particle solution
  const double s = std::max(s1, s2);
  m_Position += Direction*s;
  m_GasExit = m_Position;
  m_Location = m_RadiatorCell->IsInsideCell(m_Position) ?
               Location::Mirror : Location::MissedMirror;
}

bool ParticleTrack::SwapRadiatorCell() {
  if(m_RadiatorCell->IsInsideCell(m_Position)) {
    return true;
  }
  const auto CellNumber = m_RadiatorCell->GetCellNumber();
  if(CellNumber.first == 0 && CellNumber.second == 0) {
    return false;
  }
  auto CurrentRadiatorPosition = m_RadiatorCell->GetRadiatorPosition();
  if(m_Position.X() > m_RadiatorCell->GetHexagonSize()/2.0) {
    if(CellNumber.second > 9) {
      return false;
    }
    m_RadiatorCell++;
  } else if(m_Position.X() < -m_RadiatorCell->GetHexagonSize()/2.0) {
    if((CellNumber.first == 1 && CellNumber.second == 0) ||
       (CellNumber.first == 2 && CellNumber.second == 1)) {
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
  if(m_Location != Location::Mirror) {
    throw std::runtime_error("Cannot generate photons from tracks not at the mirror");
  }
  Photon photon = GeneratePhoton(m_AerogelEntry, m_AerogelExit, Photon::Radiator::Aerogel);
  return photon;
}

Photon ParticleTrack::GeneratePhotonFromGas() const {
  if(m_Location != Location::Mirror) {
    throw std::runtime_error("Cannot generate photons from tracks not at the mirror");
  }
  Photon photon = GeneratePhoton(m_GasEntry, m_GasExit, Photon::Radiator::Gas);
  return photon;
}

std::vector<Photon> ParticleTrack::GeneratePhotonsFromAerogel() const {
  if(m_Location != Location::Mirror) {
    throw std::runtime_error("Cannot generate photons from tracks not at the mirror");
  }
  const double RadiatorDistance = TMath::Sqrt((m_AerogelExit - m_AerogelEntry).Mag2());
  const double IndexRefraction = GetIndexRefraction(Photon::Radiator::Aerogel, 1239.8/400.0);
  const int PhotonYield = std::round(GetPhotonYield(RadiatorDistance, Beta(), IndexRefraction));
  std::vector<Photon> Photons;
  for(int i = 0; i < PhotonYield; i++) {
    Photons.push_back(GeneratePhoton(m_AerogelEntry, m_AerogelExit, Photon::Radiator::Aerogel));
  }
  return Photons;
}

std::vector<Photon> ParticleTrack::GeneratePhotonsFromGas() const {
  if(m_Location != Location::Mirror) {
    throw std::runtime_error("Cannot generate photons from tracks not at the mirror");
  }
  const double RadiatorDistance = TMath::Sqrt((m_GasExit - m_GasEntry).Mag2());
  const double IndexRefraction = GetIndexRefraction(Photon::Radiator::Gas, 1239.8/400.0);
  const int PhotonYield = std::round(GetPhotonYield(RadiatorDistance, Beta(), IndexRefraction));
  std::vector<Photon> Photons;
  Photons.reserve(PhotonYield);
  for(int i = 0; i < PhotonYield; i++) {
    Photons.emplace_back(GeneratePhoton(m_GasEntry, m_GasExit, Photon::Radiator::Gas));
  }
  return Photons;
}

double ParticleTrack::GetIndexRefraction(Photon::Radiator Radiator, double Energy) const {
  switch(Radiator) {
    case Photon::Radiator::Aerogel:
      // TODO: Add aerogel dispersion
      return 1.03;
    case Photon::Radiator::Gas:
      {
	auto GetIndex = [] (double L) {
	  return 1.0 + 0.25324*1e-6/((1.0/(73.7*73.7)) - (1.0/(L*L)));
	};
	if(m_ChromaticDispersion) {
	  // Pressure at 1.0 bar
	  // Sellmeier equation with coefficients from
	  // https://twiki.cern.ch/twiki/bin/view/LHCb/C4F10
	  // They are similar to A. Filippas, et al. Nucl. Instr. and Meth. B, 196 (2002),
	  // p. 340 but now quite...?
	  const double Lambda = 1239.841987427/Energy;
	  return GetIndex(Lambda);
	} else {
	  return GetIndex(400);
	}
      }
    default:
      return 1.0;
  }
}

Photon ParticleTrack::GeneratePhoton(const Vector &Entry,
				     const Vector &Exit,
				     Photon::Radiator Radiator) const {
  const double Energy = gRandom->Uniform(1.55, 4.31);
  const double n_phase = GetIndexRefraction(Radiator, Energy);
  const double RandomFraction = m_RandomEmissionPoint
                              ? gRandom->Uniform(0.005, 0.995) : 0.5;
  const Vector EmissionPoint = Entry + (Exit - Entry)*RandomFraction;
  const double phi = gRandom->Uniform(0.0, 2*TMath::Pi());
  const double CosTheta = 1.0/(Beta()*n_phase);
  const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
  const double CosPhi = TMath::Cos(phi);
  const double SinPhi = (phi > TMath::Pi() ? -1.0 : +1.0)
                        *TMath::Sqrt(1.0 - CosPhi*CosPhi);
  Vector Direction(SinTheta*CosPhi, SinTheta*SinPhi, CosTheta);
  // Rotate to particle frame
  const ROOT::Math::RotationY RotateY(m_Momentum.Theta());
  Direction = RotateY(Direction);
  const ROOT::Math::RotationZ RotateZ(m_Momentum.Phi());
  Direction = RotateZ(Direction);
  return {EmissionPoint,
          Direction,
          Energy,
          CosTheta,
          Radiator,
          &(*m_RadiatorCell)};
}

double ParticleTrack::Beta() const {
  const double Momentum = TMath::Sqrt(m_Momentum.Mag2());
  return Momentum/TMath::Sqrt(m_Mass*m_Mass + Momentum*Momentum);
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
  const double Efficiency = 0.60*0.90*0.80;
  const double DeltaE = 2.55;
  return x*DeltaE*37000.0*(1.0 - 1.0/TMath::Power(Beta*n, 2))*Efficiency;
}

void ParticleTrack::MapPhi(double DeltaPhi) {
  const ROOT::Math::RotationZ RotateZ(DeltaPhi);
  m_Position = RotateZ(m_Position);
  m_InitialPosition = RotateZ(m_InitialPosition);
  m_Momentum = RotateZ(m_Momentum);
}

void ParticleTrack::ReflectZ() {
  m_Position.SetZ(-m_Position.Z());
  m_InitialPosition.SetZ(-m_InitialPosition.Z());
  m_Momentum.SetZ(-m_Momentum.Z());
}

void ParticleTrack::SwapXZ(Vector &Vec) const {
  const double Temp = Vec.X();
  Vec.SetX(Vec.Z());
  Vec.SetZ(Temp);
}

std::unique_ptr<TLine> ParticleTrack::DrawParticleTrack() const {
  const auto CurrentPosition = m_RadiatorCell->GetRadiatorPosition() + m_Position;
  TLine Track(m_InitialPosition.X(),
	      m_InitialPosition.Z(),
	      CurrentPosition.X(),
	      CurrentPosition.Z());
  Track.SetLineColor(kRed);
  return std::make_unique<TLine>(Track);
}

const Vector& ParticleTrack::GetPosition() const {
  return m_Position;
}

ParticleTrack::Location ParticleTrack::GetParticleLocation() const {
  return m_Location;
}

const Vector ParticleTrack::GetEntranceWindowPosition() const {
  return m_EntranceWindowPosition;
}

double ParticleTrack::GetRadiatorColumnNumber() const {
  return m_RadiatorCell->GetCellNumber().first;
}

double ParticleTrack::GetRadiatorRowNumber() const {
  return m_RadiatorCell->GetCellNumber().second;
}
