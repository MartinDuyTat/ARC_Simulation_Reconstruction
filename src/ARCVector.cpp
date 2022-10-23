// Martin Duy Tat 16th October 2022

#include<utility>
#include<optional>
#include<stdexcept>
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Math/Rotation3D.h"
#include"Math/RotationZ.h"
#include"ARCVector.h"

using RotationZ = ROOT::Math::RotationZ;

ARCVector::ARCVector(const Vector &Vec): m_Vector(Vec) {
}

ARCVector::ARCVector(const Vector &Vec,
		     const Vector &RadiatorPosition,
		     const Rotation3D &RadiatorRotation):
  m_Vector(Vec),
  m_CoordinateTransform(std::make_pair(RadiatorPosition, RadiatorRotation)) {
}

void ARCVector::AssignLocalCoordinates(const Vector &Position,
				       const Rotation3D &Rotation) {
  if(m_CoordinateTransform.has_value()) {
    throw std::runtime_error("Cannot assign local coordinates twice");
  }
  m_CoordinateTransform = std::make_pair(Position, Rotation);
  GlobalToLocal(m_Vector);
}

bool ARCVector::IsLocal() const {
  return m_CoordinateTransform.has_value();
}

void ARCVector::ConvertToGlobal() {
  LocalToGlobal(m_Vector);
  m_CoordinateTransform.reset();
}

const Vector& ARCVector::LocalVector() const {
  if(m_CoordinateTransform.has_value()) {
    return m_Vector;
  } else {
    throw std::logic_error("Local coordinates don't exist, no local vector exists");
  }
}

Vector ARCVector::GlobalVector() const {
  if(m_CoordinateTransform.has_value()) {
    Vector Vec = m_Vector;
    LocalToGlobal(Vec);
    return Vec;
  } else {
    return m_Vector;
  }
}

ARCVector& ARCVector::operator +=(const Vector &Vec) {
  m_Vector += Vec;
  return *this;
}

void ARCVector::MapPhi(double DeltaPhi) {
  const RotationZ PhiRotation(DeltaPhi);
  if(m_CoordinateTransform.has_value()) {
    LocalToGlobal(m_Vector);
    m_Vector = PhiRotation(m_Vector);
    GlobalToLocal(m_Vector);
  } else {
    m_Vector = PhiRotation(m_Vector);
  }
}

void ARCVector::ReflectZ() {
  if(m_CoordinateTransform.has_value()) {
    LocalToGlobal(m_Vector);
    m_Vector.SetZ(-m_Vector.Z());
    GlobalToLocal(m_Vector);
  } else {
    m_Vector.SetZ(-m_Vector.Z());
  }
}

void ARCVector::ReflectY() {
  if(m_CoordinateTransform.has_value()) {
    LocalToGlobal(m_Vector);
    m_Vector.SetY(-m_Vector.Y());
    GlobalToLocal(m_Vector);
  } else {
    m_Vector.SetY(-m_Vector.Y());
  }
}

void ARCVector::SetX(double x) {
  m_Vector.SetX(x);
}

void ARCVector::SetZ(double z) {
  m_Vector.SetZ(z);
}

void ARCVector::GlobalToLocal(Vector &Vec) const {
  Vec -= m_CoordinateTransform->first;
  Vec = m_CoordinateTransform->second(Vec);
}

void ARCVector::LocalToGlobal(Vector &Vec) const {
  const auto InverseRotation = m_CoordinateTransform->second.Inverse();
  Vec = InverseRotation(Vec);
  Vec += m_CoordinateTransform->first;
}
