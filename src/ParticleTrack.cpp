// Martin Duy Tat 29th April 2022

#include<stdexcept>
#include<memory>
#include<utility>
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
#include"Utilities.h"

ParticleTrack::ParticleTrack(int ParticleID,
			     const Vector &Momentum,
			     std::size_t TrackNumber):
  Particle::Particle({}),
  m_Momentum(Momentum),
  m_InitialPosition(Vector(0.0, 0.0, 0.0)),
  m_ParticleID(ParticleID),
  m_Location(Location::TrackerVolume),
  m_RandomEmissionPoint(Settings::GetBool("General/RandomEmissionPoint")),
  m_ChromaticDispersion(Settings::GetBool("General/ChromaticDispersion")),
  m_Mass(ParticleMass::GetMass(m_ParticleID)),
  m_PhotonMultiplier(GetPhotonMultiplier(m_Momentum.GlobalVector())),
  m_TrackNumber(TrackNumber),
  m_TracksToDraw(Settings::GetSizeTVector("General/TrackToDraw")) {
}

void ParticleTrack::TrackThroughTracker(const TrackingVolume &InnerTracker) {
  // TODO: Account for magnetic field
  if(m_Location != Location::TrackerVolume) {
    throw std::runtime_error("Cannot track particle through inner tracker again");
  }
  const auto Position = m_Position.GlobalVector();
  const auto Momentum = m_Momentum.GlobalVector();
  // Check if we're considering the end cap or barrel
  if(Settings::GetString("General/BarrelOrEndcap") == "Barrel") {
    // Solve quadratic to track particle to edge of tracking barrel
    const auto MomentumUnit = Momentum.Unit();
    const double b = Position.X()*MomentumUnit.X()
                   + Position.Y()*MomentumUnit.Y();
    const double c = Position.X()*Position.X()
                   + Position.Y()*Position.Y()
                   - InnerTracker.GetRadius()*InnerTracker.GetRadius();
    const double s = TMath::Sqrt(b*b - c) - b;
    const double xyUnit = TMath::Sqrt(TMath::Power(MomentumUnit.X(), 2) +
				      TMath::Power(MomentumUnit.Y(), 2));
    m_Position += s*MomentumUnit/xyUnit;
  } else {
    // z-distance to end cap
    const double BarrelZ = Settings::GetDouble("ARCGeometry/BarrelZ");
    const double ZDist = BarrelZ - Position.Z();
    const double Slope = TMath::Abs(1.0/Momentum.Z());
    m_Position += Momentum*ZDist*Slope;
  }
  m_Location = Location::EntranceWindow;
}

void ParticleTrack::SetRadiator(const RadiatorCell *radiatorCell) {
  m_RadiatorCell = radiatorCell;
  
}

void ParticleTrack::ConvertToRadiatorCoordinates() {
  Particle::ConvertToRadiatorCoordinates();
  // Assign local coordinates
  const auto RadiatorPosition = m_RadiatorCell->GetRadiatorPosition();
  const auto RadiatorRotation = m_RadiatorCell->GetRadiatorRotation();
  m_Momentum.AssignLocalCoordinates({}, RadiatorRotation);
  m_InitialPosition.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
  m_AerogelEntry.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
  m_AerogelExit.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
  m_GasEntry.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
  m_GasExit.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
  m_EntranceWindowPosition.AssignLocalCoordinates(RadiatorPosition,
						  RadiatorRotation);
  // Check if particle is within acceptance
  if(!m_RadiatorCell->IsInsideCell(m_Position) ||
     m_Momentum.GlobalVector().Dot(m_Position.GlobalVector()) < 0.0) {
    m_Location = Location::MissedEntranceWindow;
  }
  // Save the entrance window position
  m_EntranceWindowPosition = m_Position;
}

void ParticleTrack::ConvertBackToGlobalCoordinates() {
  Particle::ConvertBackToGlobalCoordinates();
  // Remove the local coordinate system
  m_Momentum.ConvertToGlobal();
  m_InitialPosition.ConvertToGlobal();
  m_AerogelEntry.ConvertToGlobal();
  m_AerogelExit.ConvertToGlobal();
  m_GasEntry.ConvertToGlobal();
  m_GasExit.ConvertToGlobal();
  m_EntranceWindowPosition.ConvertToGlobal();
}

void ParticleTrack::TrackThroughAerogel() {
  if(!m_RadiatorCell) {
    throw std::runtime_error("No radiator cell to track through");
  }
  if(m_Location != Location::EntranceWindow) {
    throw std::runtime_error("Particle not at entrance window");
  }
  const double ZDistToAerogel = -m_Position.LocalVector().Z();
  const double Slope = 1.0/TMath::Cos(m_Momentum.LocalVector().Theta());
  const auto MomentumUnit = m_Momentum.LocalVector().Unit();
  m_Position += MomentumUnit*Slope*ZDistToAerogel;
  m_AerogelEntry = m_Position;
  m_Position += MomentumUnit*Slope*m_RadiatorCell->GetAerogelThickness();
  m_AerogelExit = m_Position;
  m_GasEntry = m_AerogelExit;
  m_Location = Location::Radiator;
}

void ParticleTrack::TrackThroughGasToMirror() {
  if(m_Location != Location::Radiator &&
     m_Location != Location::MissedEntranceWindow) {
    throw std::runtime_error("Cannot track through gas if not in radiator");
  }
  m_Position = m_GasEntry;
  const auto Position = m_GasEntry;
  const auto LocalPosition = Position.LocalVector();
  // Solve quadratic s^2 - 2sb + c = 0 to find interserction of particle track and mirror
  auto MirrorCentre = m_RadiatorCell->GetMirrorCentre();
  const double MirrorRadius = m_RadiatorCell->GetMirrorCurvature();
  const auto Direction = m_Momentum.LocalVector().Unit();
  const double b = (MirrorCentre - LocalPosition).Dot(Direction);
  const double c = MirrorCentre.Mag2()
                 + LocalPosition.Mag2()
                 - MirrorRadius*MirrorRadius
                 - 2*LocalPosition.Dot(MirrorCentre);
  const double Discriminant = b*b - c;
  if(Discriminant < 0) {
    m_Location = Location::MissedMirror;
    return;
  }
  const double s1 = b + TMath::Sqrt(Discriminant);
  const double s2 = b - TMath::Sqrt(Discriminant);
  // Pick forwards moving particle solution
  const double s = std::max(s1, s2);
  m_Position += Direction*s;
  m_GasExit = m_Position;
  if(m_RadiatorCell->IsInsideCell(m_Position)) {
    m_Location = Location::Mirror;
  }
}

bool ParticleTrack::TrackToNextCell(const RadiatorArray &radiatorArray) {
  std::size_t Counter = 0;
  while(m_Location != Location::Mirror) {
    ConvertBackToGlobalCoordinates();
    m_Location = Location::Radiator;
    if(!FindRadiator(radiatorArray) || Counter > 10) {
      return false;
    }
    TrackThroughGasToMirror();
    Counter++;
  }
  return true;
}

Photon ParticleTrack::GeneratePhotonFromAerogel() const {
  if(m_Location != Location::Radiator) {
    throw std::runtime_error("Cannot generate photons before aerogel or after gas");
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
  if(m_Location != Location::Radiator) {
    throw std::runtime_error("Cannot generate photons before aerogel or after gas");
  }
  const auto AerogelVec = m_AerogelExit.LocalVector() - m_AerogelEntry.LocalVector();
  const double RadiatorDistance = TMath::Sqrt(AerogelVec.Mag2());
  const double IndexRefraction =
    Utilities::GetIndexRefraction(Photon::Radiator::Aerogel, false);
  const std::size_t PhotonYield =static_cast<std::size_t>(
    std::round(GetPhotonYield(RadiatorDistance, Beta(), IndexRefraction)));
  std::vector<Photon> Photons;
  for(std::size_t i = 0; i < PhotonYield; i++) {
    Photons.push_back(GeneratePhoton(m_AerogelEntry,
				     m_AerogelExit,
				     Photon::Radiator::Aerogel));
  }
  return Photons;
}

std::vector<Photon> ParticleTrack::GeneratePhotonsFromGas() const {
  if(m_Location != Location::Mirror) {
    throw std::runtime_error("Cannot generate photons from tracks not at the mirror");
  }
  const auto GasVec = m_GasExit.LocalVector() - m_GasEntry.LocalVector();
  const double RadiatorDistance = TMath::Sqrt(GasVec.Mag2());
  const double IndexRefraction =
    Utilities::GetIndexRefraction(Photon::Radiator::Gas, false);
  const std::size_t PhotonYield = static_cast<std::size_t>(
    std::round(GetPhotonYield(RadiatorDistance, Beta(), IndexRefraction)));
  std::vector<Photon> Photons;
  Photons.reserve(PhotonYield);
  for(std::size_t i = 0; i < PhotonYield; i++) {
    Photons.emplace_back(GeneratePhoton(m_GasEntry, m_GasExit, Photon::Radiator::Gas));
  }
  return Photons;
}

Photon ParticleTrack::GeneratePhoton(const ARCVector &Entry,
				     const ARCVector &Exit,
				     Photon::Radiator Radiator) const {
  const auto EntryVec = Entry.LocalVector();
  const auto ExitVec = Exit.LocalVector();
  const double Energy = gRandom->Uniform(1.55, 4.31);
  const double n_phase = Utilities::GetIndexRefraction(Radiator,
						       m_ChromaticDispersion,
						       Energy);
  const double RandomFraction = m_RandomEmissionPoint
                              ? gRandom->Uniform(0.005, 0.995) : 0.5;
  const Vector EmissionPoint = EntryVec + (ExitVec - EntryVec)*RandomFraction;
  const Vector AssumedEmissionPoint = 0.5*(EntryVec + ExitVec);
  const double phi = gRandom->Uniform(0.0, 2*TMath::Pi());
  auto GetCosTheta = [&] () {
    double CosTheta = 1.0/(Beta()*n_phase);
    if(CosTheta > 1.0) {
      return 1.0;
    } else {
      return CosTheta;
    }
  };
  const double CosTheta = GetCosTheta();
  const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
  const double CosPhi = TMath::Cos(phi);
  const double SinPhi = (phi > TMath::Pi() ? -1.0 : +1.0)
                        *TMath::Sqrt(1.0 - CosPhi*CosPhi);
  Vector Direction(SinTheta*CosPhi, SinTheta*SinPhi, CosTheta);
  // Rotate to particle frame
  const ROOT::Math::RotationY RotateY(m_Momentum.LocalVector().Theta());
  Direction = RotateY(Direction);
  const ROOT::Math::RotationZ RotateZ(m_Momentum.LocalVector().Phi());
  Direction = RotateZ(Direction);
  const auto RadiatorPosition = m_RadiatorCell->GetRadiatorPosition();
  const auto RadiatorRotation = m_RadiatorCell->GetRadiatorRotation();
  Photon photon(ARCVector(EmissionPoint, RadiatorPosition, RadiatorRotation),
		ARCVector(AssumedEmissionPoint, RadiatorPosition, RadiatorRotation),
		ARCVector(Direction, {}, RadiatorRotation),
		m_Momentum.GlobalVector().Unit(),
		Energy,
		CosTheta,
		Radiator,
		&(*m_RadiatorCell),
		IsTrackDrawn(),
		1.0/m_PhotonMultiplier);
  return photon;
}

double ParticleTrack::Beta() const {
  const double Momentum = TMath::Sqrt(m_Momentum.LocalVector().Mag2());
  return Momentum/TMath::Sqrt(m_Mass*m_Mass + Momentum*Momentum);
}

const Vector& ParticleTrack::GetEntryPoint(Photon::Radiator Radiator) const {
  if(Radiator == Photon::Radiator::Gas) {
    return m_GasEntry.LocalVector();
  } else if(Radiator == Photon::Radiator::Aerogel) {
    return m_AerogelEntry.LocalVector();
  } else {
    throw std::runtime_error("Cannot find entry point of unknown radiator");
  }
}

const Vector& ParticleTrack::GetExitPoint(Photon::Radiator Radiator) const {
  if(Radiator == Photon::Radiator::Gas) {
    return m_GasExit.LocalVector();
  } else if(Radiator == Photon::Radiator::Aerogel) {
    return m_AerogelExit.LocalVector();
  } else {
    throw std::runtime_error("Cannot find entry point of unknown radiator");
  }
}

double ParticleTrack::GetPhotonYield(double x, double Beta, double n) const {
  // TODO: Move this to separate class
  if(Beta*n < 1.0) {
    return 0.0;
  } else {
    const double Efficiency = 0.60*0.90*0.80;
    const double DeltaE = 2.55;
    const double SinThetaC2 = (1.0 - 1.0/TMath::Power(Beta*n, 2));
    return x*DeltaE*37000.0*SinThetaC2*Efficiency*m_PhotonMultiplier;
  }
}

void ParticleTrack::MapPhi(double DeltaPhi) {
  Particle::MapPhi(DeltaPhi);
  m_InitialPosition.MapPhi(DeltaPhi);
  m_Momentum.MapPhi(DeltaPhi);
  m_AerogelEntry.MapPhi(DeltaPhi);
  m_AerogelExit.MapPhi(DeltaPhi);
  m_GasEntry.MapPhi(DeltaPhi);
  m_GasExit.MapPhi(DeltaPhi);
}

void ParticleTrack::ReflectZ() {
  Particle::ReflectZ();
  m_InitialPosition.ReflectZ();
  m_Momentum.ReflectZ();
  m_AerogelEntry.ReflectZ();
  m_AerogelExit.ReflectZ();
  m_GasEntry.ReflectZ();
  m_GasExit.ReflectZ();
}

void ParticleTrack::ReflectY() {
  Particle::ReflectY();
  m_InitialPosition.ReflectY();
  m_Momentum.ReflectY();
  m_AerogelEntry.ReflectY();
  m_AerogelExit.ReflectY();
  m_GasEntry.ReflectY();
  m_GasExit.ReflectY();
}

std::unique_ptr<TLine> ParticleTrack::DrawParticleTrack() const {
  if(IsTrackDrawn()) {
    const auto ReversePhi = m_RadiatorCell->ReversePhiRotation();
    const auto InitialPosition = ReversePhi(m_InitialPosition.GlobalVector());
    const auto Position = ReversePhi(m_Position.GlobalVector());
    TLine Track(InitialPosition.Z(),
		InitialPosition.X(),
		Position.Z(),
		Position.X());
    Track.SetLineColor(kRed);
    return std::make_unique<TLine>(Track);
  } else {
    return nullptr;
  }
}

ParticleTrack::Location ParticleTrack::GetParticleLocation() const {
  return m_Location;
}

Vector ParticleTrack::GetEntranceWindowPosition() const {
  return m_EntranceWindowPosition.GlobalVector();
}

bool ParticleTrack::IsAtRadiator() const {
  if(m_Location != Location::EntranceWindow &&
     m_Location != Location::MissedEntranceWindow &&
     m_Location != Location::Radiator) {
    return false;
  } else {
    return true;
  }
}
bool ParticleTrack::IsTrackDrawn() const {
  if(!Settings::Exists("EventDisplay/RowToDraw")) {
    return false;
  }
  const std::size_t RowToDraw = Settings::GetSizeT("EventDisplay/RowToDraw");
  if(RowToDraw == GetRadiatorRowNumber()) {
    const bool DrawAllTracks = Settings::GetBool("General/DrawAllTracks");
    if(DrawAllTracks) {
      return true;
    }
    const auto iter = std::find(m_TracksToDraw.begin(),
				m_TracksToDraw.end(),
				m_TrackNumber);
    if(iter != m_TracksToDraw.end()) {
      return true;
    }
  }
  return false;
}

double ParticleTrack::GetPhotonMultiplier(const Vector &Momentum) {
  double LowMomentumLimit = Settings::GetDouble("General/LowMomentumLimit");
  double PhotonMultiplier = Settings::GetDouble("General/PhotonMultiplier");
  const std::string GasOrAerogel = Settings::GetString("General/GasOrAerogel");
  if(GasOrAerogel == "Gas" && TMath::Sqrt(Momentum.Mag2()) < LowMomentumLimit) {
    return PhotonMultiplier*20.0;
  } else {
    return PhotonMultiplier;
  }
}
