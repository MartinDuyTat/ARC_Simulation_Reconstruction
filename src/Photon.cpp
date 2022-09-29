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

Photon::Photon(const Vector &Position,
	       const Vector &Direction,
	       double Energy,
	       double CosCherenkovAngle,
	       Radiator radiator,
	       const RadiatorCell *radiatorCell):
  Particle(Position, Particle::CoordinateSystem::LocalRadiator, radiatorCell),
  m_EmissionPoint(Position),
  m_Direction(Direction),
  m_Energy(Energy),
  m_Radiator(radiator),
  m_CosCherenkovAngle(CosCherenkovAngle),
  m_Status(Status::Emitted),
  m_AerogelTravelDistance(0.0),
  m_MirrorHitPosition(nullptr) {
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
Photon::DrawPhotonPath() const {
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> PhotonLine;
  const auto RadiatorPosition = m_RadiatorCell->GetRadiatorPosition();
  const auto EmissionPoint = Utilities::SwapXZForEndCap(m_RadiatorCell,
							m_EmissionPoint
							+ RadiatorPosition);
  const auto Position = Utilities::SwapXZForEndCap(m_RadiatorCell,
						   m_Position + RadiatorPosition);
  if(!m_MirrorHitPosition) {
    TLine PhotonLine1(EmissionPoint.X(),
		      EmissionPoint.Z(),
		      Position.X(),
		      Position.Z());
    if(m_Status == Status::MirrorMiss) {
      PhotonLine1.SetLineColor(7);
    } else {
      PhotonLine1.SetLineColor(kGreen);
    }
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine1), ""));
  } else {
    const auto MirrorHitPosition = Utilities::SwapXZForEndCap(m_RadiatorCell,
							      *m_MirrorHitPosition 
							      + RadiatorPosition);
    TLine PhotonLine1(EmissionPoint.X(),
		      EmissionPoint.Z(),
		      MirrorHitPosition.X(),
		      MirrorHitPosition.Z());
    TLine PhotonLine2(MirrorHitPosition.X(),
		      MirrorHitPosition.Z(),
		      Position.X(),
		      Position.Z());
    if(m_Status != Status::DetectorHit) {
      PhotonLine1.SetLineColor(6);
      PhotonLine2.SetLineColor(6);
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
  return m_Direction;
}

void Photon::KickPhoton(const Vector &Kick) {
  m_Direction += Kick;
}

Photon::Status Photon::GetStatus() const {
  return m_Status;
}

const Vector& Photon::GetEmissionPoint() const {
  return m_EmissionPoint;
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

const Vector* Photon::GetMirrorHitPosition() const {
  return m_MirrorHitPosition.get();
}

void Photon::RegisterMirrorHitPosition(const Vector &MirrorHitPosition) {
  if(m_MirrorHitPosition) {
    throw std::runtime_error("Mirror hit position already exists");
  } else {
    m_MirrorHitPosition = std::make_unique<Vector>(MirrorHitPosition);
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
  // Check if particle is within acceptance
  if(!m_RadiatorCell->IsInsideCell(m_Position)) {
    m_Status = Status::Outside;
  }
}

void Photon::ConvertBackToGlobalCoordinates() {
  Particle::ConvertBackToGlobalCoordinates();
  // Photon doesn't hit any mirror
  m_MirrorHitPosition = nullptr;
}

void Photon::MapPhi(double DeltaPhi) {
  Particle::MapPhi(DeltaPhi);
  const ROOT::Math::RotationZ RotateZ(DeltaPhi);
  m_EmissionPoint = RotateZ(m_EmissionPoint);
  m_Direction = RotateZ(m_Direction);
}

void Photon::ReflectZ() {
  Particle::ReflectZ();
  m_EmissionPoint.SetZ(-m_EmissionPoint.Z());
  m_Direction.SetZ(m_Direction.Z());
}

void Photon::ReflectY() {
  Particle::ReflectY();
  m_EmissionPoint.SetY(-m_EmissionPoint.Y());
  m_Direction.SetY(m_Direction.Y());
}

void Photon::SwapXZ() {
  Particle::SwapXZ();
  Particle::SwapXZ(m_EmissionPoint);
  Particle::SwapXZ(m_Direction);
}

bool Photon::IsAtRadiator() const {
  return m_Status == Status::MirrorMiss;
}

void Photon::ChangeCoordinateOrigin(const Vector &Shift) {
  Particle::ChangeCoordinateOrigin(Shift);
  m_EmissionPoint -= Shift;
}
