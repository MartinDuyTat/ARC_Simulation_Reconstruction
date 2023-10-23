// Martin Duy Tat 29th April 2022

#include<memory>
#include<vector>
#include<utility>
#include<string>
#include<stdexcept>
#include"TLine.h"
#include"Math/RotationY.h"
#include"Math/RotationZ.h"
#include"Photon.h"
#include"RadiatorCell.h"
#include"BarrelRadiatorCell.h"
#include"Utilities.h"
#include"ARCVector.h"
#include"Settings.h"

Photon::Photon(const Vector &Position,
	       const Vector &AssumedPosition,
	       const Vector &Direction,
	       const Vector &ParticleDirection,
	       const Vector &AssumedParticleDirection,
	       double Energy,
	       double CosCherenkovAngle,
	       Radiator radiator,
	       const RadiatorCell *radiatorCell, 
	       bool IsTrackDrawn,
	       double Weight):
  Particle(Position, radiatorCell),
  m_AssumedEmissionPoint(AssumedPosition),
  m_EmissionPoint(Position),
  m_Direction(Direction),
  m_ParticleDirection(ParticleDirection),
  m_AssumedParticleDirection(AssumedParticleDirection),
  m_Energy(Energy),
  m_Radiator(radiator),
  m_CosCherenkovAngle(CosCherenkovAngle),
  m_Status(Status::Emitted),
  m_AerogelTravelDistance(0.0),
  m_MirrorHitPosition(nullptr),
  m_Weight(Weight),
  m_HasMigrated(false),
  m_IsTrackDrawn(IsTrackDrawn) {
  const auto RadiatorPosition = radiatorCell->GetRadiatorPosition();
  const auto RadiatorRotation = radiatorCell->GetRadiatorRotation();
  m_Position.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
  m_AssumedEmissionPoint.AssignLocalCoordinates(RadiatorPosition,
						RadiatorRotation);
  m_EmissionPoint.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
  m_Direction.AssignLocalCoordinates({}, RadiatorRotation);
  m_ParticleDirection.AssignLocalCoordinates({}, RadiatorRotation);
  m_AssumedParticleDirection.AssignLocalCoordinates({}, RadiatorRotation);
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
Photon::DrawPhotonPath() const {
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> PhotonLine;
  if(!m_IsTrackDrawn) {
    return PhotonLine;
  }
  const std::size_t RowToDraw = Settings::GetSizeT("EventDisplay/RowToDraw");
  if(!m_RadiatorCell || GetRadiatorRowNumber() != RowToDraw) {
    return PhotonLine;
  }
  const auto ReversePhi = m_RadiatorCell->ReversePhiRotation();
  const auto EmissionPoint = ReversePhi(m_EmissionPoint.GlobalVector());
  const auto Position = ReversePhi(m_Position.GlobalVector());
  if(!m_MirrorHitPosition) {
    TLine PhotonLine1(EmissionPoint.Z(),
		      EmissionPoint.X(),
		      Position.Z(),
		      Position.X());
    if(m_Status == Status::MirrorMiss) {
      PhotonLine1.SetLineColor(7);
    } else {
      PhotonLine1.SetLineColor(kGreen);
    }
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine1), ""));
  } else {
    const auto MirrorHitPosition = ReversePhi(m_MirrorHitPosition->GlobalVector());
    TLine PhotonLine1(EmissionPoint.Z(),
		      EmissionPoint.X(),
		      MirrorHitPosition.Z(),
		      MirrorHitPosition.X());
    TLine PhotonLine2(MirrorHitPosition.Z(),
		      MirrorHitPosition.X(),
		      Position.Z(),
		      Position.X());
    if(m_Status != Status::DetectorHit) {
      PhotonLine1.SetLineColor(6);
      PhotonLine2.SetLineColor(6);
    } else if(m_Status == Status::DetectorMiss) {
      PhotonLine1.SetLineColor(8);
      PhotonLine2.SetLineColor(8);
    } else {
      PhotonLine1.SetLineColor(kBlue);
      PhotonLine2.SetLineColor(kBlue);
    }
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine1), ""));
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine2), ""));
  }
  return PhotonLine;
}

void Photon::PropagatePhoton(const Vector &Displacement) {
  m_Position += Displacement;
}

const Vector& Photon::GetDirection() const {
  return m_Direction.LocalVector();
}

void Photon::KickPhoton(const Vector &Kick) {
  m_Direction += Kick;
}

Photon::Status Photon::GetStatus() const {
  return m_Status;
}

const Vector& Photon::GetEmissionPoint(bool TrueEmissionPoint) const {
  if(TrueEmissionPoint) {
    return m_EmissionPoint.LocalVector();
  } else {
    return m_AssumedEmissionPoint.LocalVector();
  }
}

void Photon::UpdatePhotonStatus(Photon::Status status) {
  if(m_Status == status) {
    throw std::runtime_error("New photon status must be different");
  }
  m_Status = status;
}

double Photon::GetEnergy() const {
  return m_Energy;
}

void Photon::AddAerogelTravelDistance(double Distance) {
  if(Distance < 0.0) {
    throw std::runtime_error("Aerogel travel distance cannot be negative");
  }
  m_AerogelTravelDistance += Distance;
}

double Photon::GetAerogelTravelDistance() const {
  return m_AerogelTravelDistance;
}

const ARCVector* Photon::GetMirrorHitPosition() const {
  return m_MirrorHitPosition.get();
}

void Photon::RegisterMirrorHitPosition(const Vector &MirrorHitPosition) {
  if(m_MirrorHitPosition) {
    throw std::runtime_error("Mirror hit position already exists");
  } else {
    m_MirrorHitPosition = std::make_unique<ARCVector>(MirrorHitPosition);
  }
}

Photon::Radiator Photon::GetRadiator() const {
  return m_Radiator;
}

double Photon::GetCosCherenkovAngle() const {
  return m_CosCherenkovAngle;
}

void Photon::ConvertToRadiatorCoordinates() {
  Particle::ConvertToRadiatorCoordinates();
  const auto RadiatorPosition = m_RadiatorCell->GetRadiatorPosition();
  const auto RadiatorRotation = m_RadiatorCell->GetRadiatorRotation();
  m_AssumedEmissionPoint.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
  m_EmissionPoint.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
  m_Direction.AssignLocalCoordinates({}, RadiatorRotation);
  m_ParticleDirection.AssignLocalCoordinates({}, RadiatorRotation);
  m_AssumedParticleDirection.AssignLocalCoordinates({}, RadiatorRotation);
  // Check if particle is within acceptance
  if(!m_RadiatorCell->IsInsideCell(m_Position)) {
    m_Status = Status::OutsideCell;
  }
}

void Photon::ConvertBackToGlobalCoordinates() {
  Particle::ConvertBackToGlobalCoordinates();
  m_AssumedEmissionPoint.ConvertToGlobal();
  m_EmissionPoint.ConvertToGlobal();
  m_Direction.ConvertToGlobal();
  m_ParticleDirection.ConvertToGlobal();
  m_AssumedParticleDirection.ConvertToGlobal();
  // Photon doesn't hit any mirror
  m_MirrorHitPosition = nullptr;
}

void Photon::MapPhi(double DeltaPhi) {
  Particle::MapPhi(DeltaPhi);
  m_AssumedEmissionPoint.MapPhi(DeltaPhi);
  m_EmissionPoint.MapPhi(DeltaPhi);
  m_Direction.MapPhi(DeltaPhi);
  m_ParticleDirection.MapPhi(DeltaPhi);
  m_AssumedParticleDirection.MapPhi(DeltaPhi);
}

void Photon::ReflectZ() {
  Particle::ReflectZ();
  m_AssumedEmissionPoint.ReflectZ();
  m_EmissionPoint.ReflectZ();
  m_Direction.ReflectZ();
  m_ParticleDirection.ReflectZ();
  m_AssumedParticleDirection.ReflectZ();
}


void Photon::ReflectY() {
  Particle::ReflectY();
  m_AssumedEmissionPoint.ReflectY();
  m_EmissionPoint.ReflectY();
  m_Direction.ReflectY();;
  m_ParticleDirection.ReflectY();
  m_AssumedParticleDirection.ReflectY();
}

bool Photon::IsAtRadiator() const {
  return m_Status == Status::Emitted;
}

void Photon::PutPhotonToEmissionPoint() {
  m_Position = m_EmissionPoint;
}

void Photon::PhotonHasMigrated() {
  m_HasMigrated = true;
}

bool Photon::HasPhotonMigrated() const {
  return m_HasMigrated;
}

const Vector& Photon::GetParticleDirection(bool TrueEmissionPoint) const {
  if(TrueEmissionPoint) {
    return m_ParticleDirection.LocalVector();
  } else {
    return m_AssumedParticleDirection.LocalVector();
  }
}

double Photon::GetWeight() const {
  return m_Weight;
}
