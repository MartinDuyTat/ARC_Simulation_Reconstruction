// Martin Duy Tat 29th September 2022

#include"Math/RotationY.h"
#include"Math/RotationZ.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Particle.h"
#include"RadiatorCell.h"
#include"Settings.h"

Particle::Particle(const Vector &Position,
		   CoordinateSystem coordinateSystem,
		   const RadiatorCell *radiatorCell):
  m_Position(Position),
  m_RadiatorCell(radiatorCell),
  m_CoordinateSystem(coordinateSystem),
  m_PhiRotated(0.0) {
}

const Vector& Particle::GetPosition() const {
  return m_Position;
}


void Particle::SetPhiRotated(double Phi) {
  m_PhiRotated = Phi;
}

const RadiatorCell* Particle::GetRadiatorCell() const {
  return m_RadiatorCell;
}

void Particle::SwapXZ(Vector &Vec) const {
  const double Temp = Vec.X();
  Vec.SetX(Vec.Z());
  Vec.SetZ(Temp);
}
void Particle::ConvertToRadiatorCoordinates() {
  // Check if coordinate system if correct
  if(m_CoordinateSystem == CoordinateSystem::LocalRadiator) {
    throw std::runtime_error("Particle position is already in local radiator coordinates");
  }
  // For barrel, rotate around y-axis so that z axis now points towards the high pT cell
  if(Settings::GetString("General/BarrelOrEndcap") == "Barrel") {
    SwapXZ();
  }
  ChangeCoordinateOrigin(m_RadiatorCell->GetRadiatorPosition());
  m_CoordinateSystem = CoordinateSystem::LocalRadiator;
}

void Particle::ConvertBackToGlobalCoordinates() {
  // Check if coordinate system if correct
  if(m_CoordinateSystem == CoordinateSystem::GlobalDetector) {
    throw std::runtime_error("Particle position is already in global detector coordinates");
  }
  // Shift coordinates so that the origin is at the global detector origin
  ChangeCoordinateOrigin(-m_RadiatorCell->GetRadiatorPosition());
  // For barrel, rotate around y-axis so that z axis now points along the barrel
  if(Settings::GetString("General/BarrelOrEndcap") == "Barrel") {
    SwapXZ();
  }
  m_CoordinateSystem = CoordinateSystem::GlobalDetector;
  // Remove the radiator cell
  m_RadiatorCell = nullptr;
  // Rotate back in phi (for upper row cells)
  MapPhi(-m_PhiRotated);
  m_PhiRotated = 0.0;
}

void Particle::MapPhi(double DeltaPhi) {
  // Check if coordinate system if correct
  if(m_CoordinateSystem == CoordinateSystem::LocalRadiator) {
    throw std::runtime_error("Cannot rotate in phi with local coordinates");
  }
  const ROOT::Math::RotationZ RotateZ(DeltaPhi);
  m_Position = RotateZ(m_Position);
}

void Particle::ReflectZ() {
  m_Position.SetZ(-m_Position.Z());
}

void Particle::ReflectY() {
  m_Position.SetY(-m_Position.Y());
}

void Particle::SwapXZ() {
  SwapXZ(m_Position);
}

void Particle::ChangeCoordinateOrigin(const Vector &Shift) {
  m_Position -= Shift;
}
