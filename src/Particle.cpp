// Martin Duy Tat 29th September 2022

#include<stdexcept>
#include"Math/Rotation3D.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Particle.h"
#include"RadiatorCell.h"
#include"RadiatorArray.h"
#include"Settings.h"

Particle::Particle(const Vector &Position,
		   const RadiatorCell *radiatorCell):
  m_Position(Position),
  m_RadiatorCell(radiatorCell) {
}

const ARCVector& Particle::GetPosition() const {
  return m_Position;
}

const RadiatorCell* Particle::GetRadiatorCell() const {
  return m_RadiatorCell;
}

bool Particle::FindRadiator(const RadiatorArray &radiatorArray) {
  m_RadiatorCell = radiatorArray.FindRadiator(*this);
  if(m_RadiatorCell) {
    ConvertToRadiatorCoordinates();
    return true;
  } else {
    return false;
  }
}

void Particle::ConvertToRadiatorCoordinates() {
  // Check if radiator cell exists
  if(!m_RadiatorCell) {
    throw std::runtime_error("Cannot switch to non-existent local coordinates");
  }
  // Assign the local coordinate system
  const auto RadiatorPosition = m_RadiatorCell->GetRadiatorPosition();
  const auto RadiatorRotation = m_RadiatorCell->GetRadiatorRotation();
  m_Position.AssignLocalCoordinates(RadiatorPosition, RadiatorRotation);
}

void Particle::ConvertBackToGlobalCoordinates() {
  // Remove the local coordinate system
  m_Position.ConvertToGlobal();
  // Remove the radiator cell
  m_RadiatorCell = nullptr;
}

void Particle::MapPhi(double DeltaPhi) {
  m_Position.MapPhi(DeltaPhi);
}

void Particle::ReflectZ() {
  m_Position.ReflectZ();
}

void Particle::ReflectY() {
  m_Position.ReflectY();
}

std::size_t Particle::GetRadiatorColumnNumber() const {
  if(m_RadiatorCell) {
    return m_RadiatorCell->GetCellNumber().first;
  } else {
    return static_cast<std::size_t>(-1);;
  }
}

std::size_t Particle::GetRadiatorRowNumber() const {
  if(m_RadiatorCell) {
    return m_RadiatorCell->GetCellNumber().second;
  } else {
    return static_cast<std::size_t>(-1);;
  }
}
