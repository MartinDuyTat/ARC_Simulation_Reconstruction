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
#include"HelixFunctor.h"

ParticleTrack::ParticleTrack(int ParticleID,
			     const Vector &Momentum,
			     std::size_t TrackNumber,
			     double BField):
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
  m_TracksToDraw(Settings::GetSizeTVector("General/TrackToDraw")),
  m_Helix(Momentum, m_ParticleID > 0 ? +1 : -1, BField),
  m_PathLength(0.0) {
}

bool ParticleTrack::TrackThroughTracker(const TrackingVolume &InnerTracker) {
  if(m_Location != Location::TrackerVolume) {
    throw std::runtime_error("Cannot track particle through inner tracker again");
  }
  if(Settings::GetString("General/BarrelOrEndcap") == "Barrel") {
    const double TrackerRadius = InnerTracker.GetRadius();
    BarrelHelixFunctor Functor(TrackerRadius);
    const double Tracker_s = m_Helix.SolvePathLength(Functor,
						     0.9*TrackerRadius,
						     5.0*TrackerRadius);
    if(Tracker_s == -999.0) {
      // Most likely the track didn't reach the barrel
      return false;
    }
    const auto NewPosition = m_Helix.GetPosition(Tracker_s);
    m_Position.SetGlobalVector(NewPosition);
    m_PathLength = Tracker_s;
  } else {
    throw std::runtime_error("Haven't implemented end cap magnetic field tracking");
  }
  m_Location = Location::EntranceWindow;
  return true;
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
  m_EntranceWindowPosition.ConvertToGlobal();
}

bool ParticleTrack::TrackThroughAerogel() {
  if(!m_RadiatorCell) {
    throw std::runtime_error("No radiator cell to track through");
  }
  if(m_Location != Location::EntranceWindow) {
    throw std::runtime_error("Particle not at entrance window");
  }
  ZPlaneHelixFunctor CoolingFunctor(0.0,
				    m_RadiatorCell->GetRadiatorPosition(),
				    true);
  const double Cooling_s = m_Helix.SolvePathLength(CoolingFunctor,
						   m_PathLength,
						   1.2*m_PathLength);
  if(Cooling_s == -999.0) {
    return false;
  }
  const auto CoolingLayerPosition = m_Helix.GetPosition(Cooling_s);
  m_Position.SetGlobalVector(CoolingLayerPosition);
  m_PathLength = Cooling_s;
  m_AerogelEntry_s = Cooling_s;
  if(Settings::GetString("General/BarrelOrEndcap") == "Barrel") {
    //const auto Position = m_Position.GlobalVector();
    const double AerogelThickness = m_RadiatorCell->GetAerogelThickness();
    ZPlaneHelixFunctor AerogelFunctor(AerogelThickness,
				      m_RadiatorCell->GetRadiatorPosition(),
				      true);
    const double Aerogel_s = m_Helix.SolvePathLength(AerogelFunctor,
						     m_PathLength,
						     1.2*m_PathLength);
    if(Aerogel_s == -999.0) {
      return false;
    }
    const auto AerogelPosition = m_Helix.GetPosition(Aerogel_s);
    m_Position.SetGlobalVector(AerogelPosition);
    m_PathLength = Aerogel_s;
  } else {
    throw std::runtime_error("Haven't implemented end cap magnetic field tracking");
  }
  m_AerogelExit_s = m_PathLength;
  m_GasEntry_s = m_AerogelExit_s;
  m_Location = Location::Radiator;
  return true;
}

void ParticleTrack::TrackThroughGasToMirror() {
  if(m_Location != Location::Radiator &&
     m_Location != Location::MissedEntranceWindow) {
    throw std::runtime_error("Cannot track through gas if not in radiator");
  }
  m_Position.SetGlobalVector(m_Helix.GetPosition(m_GasEntry_s));
  const auto MirrorCentre = m_RadiatorCell->GetGlobalMirrorCentre();
  MirrorHelixFunctor MirrorFunctor(MirrorCentre,
				   m_RadiatorCell->GetMirrorCurvature());
  const double Mirror_s = m_Helix.SolvePathLength(MirrorFunctor,
						  0.9*m_PathLength,
						  2.0*m_PathLength);
  const auto MirrorPosition = m_Helix.GetPosition(Mirror_s);
  m_Position.SetGlobalVector(MirrorPosition);
  m_PathLength = Mirror_s;
  m_GasExit_s = Mirror_s;
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
    throw std::runtime_error(
      "Cannot generate photons before aerogel or after gas");
  }
  Photon photon = GeneratePhoton(m_AerogelEntry_s,
				 m_AerogelExit_s,
				 Photon::Radiator::Aerogel);
  return photon;
}

Photon ParticleTrack::GeneratePhotonFromGas() const {
  if(m_Location != Location::Mirror) {
    throw std::runtime_error(
      "Cannot generate photons from tracks not at the mirror");
  }
  Photon photon = GeneratePhoton(m_GasEntry_s,
				 m_GasExit_s,
				 Photon::Radiator::Gas);
  return photon;
}

std::vector<Photon> ParticleTrack::GeneratePhotonsFromAerogel() const {
  if(m_Location != Location::Radiator) {
    throw std::runtime_error(
      "Cannot generate photons before aerogel or after gas");
  }
  const double RadiatorDistance = m_AerogelExit_s - m_AerogelEntry_s;
  const double IndexRefraction =
    Utilities::GetIndexRefraction(Photon::Radiator::Aerogel, false);
  const std::size_t PhotonYield =static_cast<std::size_t>(
    std::round(GetPhotonYield(RadiatorDistance, Beta(), IndexRefraction)));
  std::vector<Photon> Photons;
  for(std::size_t i = 0; i < PhotonYield; i++) {
    Photons.push_back(GeneratePhoton(m_AerogelEntry_s,
				     m_AerogelExit_s,
				     Photon::Radiator::Aerogel));
  }
  return Photons;
}

std::vector<Photon> ParticleTrack::GeneratePhotonsFromGas() const {
  if(m_Location != Location::Mirror) {
    throw std::runtime_error(
      "Cannot generate photons from tracks not at the mirror");
  }
  const double RadiatorDistance = m_GasExit_s - m_GasEntry_s;
  const double IndexRefraction =
    Utilities::GetIndexRefraction(Photon::Radiator::Gas, false);
  const std::size_t PhotonYield = static_cast<std::size_t>(
    std::round(GetPhotonYield(RadiatorDistance, Beta(), IndexRefraction)));
  std::vector<Photon> Photons;
  Photons.reserve(PhotonYield);
  for(std::size_t i = 0; i < PhotonYield; i++) {
    Photons.emplace_back(GeneratePhoton(m_GasEntry_s,
					m_GasExit_s,
					Photon::Radiator::Gas));
  }
  return Photons;
}

Photon ParticleTrack::GeneratePhoton(double Entry_s,
				     double Exit_s,
				     Photon::Radiator Radiator) const {
  // Random emission point
  const double RandomFraction = m_RandomEmissionPoint
                              ? gRandom->Uniform(0.005, 0.995) : 0.5;
  const double EmissionPoint_s = Entry_s + (Exit_s - Entry_s)*RandomFraction;
  const auto EmissionPoint = m_Helix.GetPosition(EmissionPoint_s);
  // Assumed emission point
  const double AssumedEmissionPoint_s = 0.5*(Entry_s + Exit_s);
  const auto AssumedEmissionPoint = m_Helix.GetPosition(AssumedEmissionPoint_s);
  // Energy and index of refraction
  const double Energy = gRandom->Uniform(1.55, 4.31);
  const double n_phase = Utilities::GetIndexRefraction(Radiator,
						       m_ChromaticDispersion,
						       Energy);
  // Random azimuthal angle
  const double phi = gRandom->Uniform(0.0, 2*TMath::Pi());
  // Cherenkov angle
  auto GetCosTheta = [&] () {
    double CosTheta = 1.0/(Beta()*n_phase);
    if(CosTheta > 1.0) {
      return 1.0;
    } else {
      return CosTheta;
    }
  };
  // Trigonometry calculations
  const double CosTheta = GetCosTheta();
  const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
  const double CosPhi = TMath::Cos(phi);
  const double SinPhi = (phi > TMath::Pi() ? -1.0 : +1.0)
                        *TMath::Sqrt(1.0 - CosPhi*CosPhi);
  // Photon direction in particle rest frame
  Vector Direction(SinTheta*CosPhi, SinTheta*SinPhi, CosTheta);
  // Rotate to lab frame
  const auto ParticleDirection = m_Helix.GetDirection(EmissionPoint_s);
  const ROOT::Math::RotationY RotateY(ParticleDirection.Theta());
  Direction = RotateY(Direction);
  const ROOT::Math::RotationZ RotateZ(ParticleDirection.Phi());
  Direction = RotateZ(Direction);
  // Finally create the photon
  return Photon(EmissionPoint,
		AssumedEmissionPoint,
		Direction,
		ParticleDirection,
		m_Helix.GetDirection(AssumedEmissionPoint_s),
		Energy,
		CosTheta,
		Radiator,
		&(*m_RadiatorCell),
		IsTrackDrawn(),
		1.0/m_PhotonMultiplier);
}

double ParticleTrack::Beta() const {
  const double Momentum = TMath::Sqrt(m_Momentum.LocalVector().Mag2());
  return Momentum/TMath::Sqrt(m_Mass*m_Mass + Momentum*Momentum);
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
  m_Helix.ResetOrigin(m_PathLength);
  m_Helix.MapPhi(DeltaPhi);
}

void ParticleTrack::ReflectZ() {
  Particle::ReflectZ();
  m_InitialPosition.ReflectZ();
  m_Momentum.ReflectZ();
  m_Helix.ReflectZ();
}

void ParticleTrack::ReflectY() {
  Particle::ReflectY();
  m_InitialPosition.ReflectY();
  m_Momentum.ReflectY();
  m_Helix.ResetOrigin(m_PathLength);
  m_Helix.ReflectY();
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
