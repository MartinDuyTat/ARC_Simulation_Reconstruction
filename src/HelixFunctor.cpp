// Martin Duy Tat 28th October 2022

#include"TMath.h"
#include"HelixFunctor.h"

BarrelHelixFunctor::BarrelHelixFunctor(double Radius): m_Radius(Radius) {
}

double BarrelHelixFunctor::operator()(const Vector &Vect) const {
  return Vect.X()*Vect.X() + Vect.Y()*Vect.Y() - m_Radius*m_Radius;
}

EndCapHelixFunctor::EndCapHelixFunctor(double z): m_z(z) {
}

double EndCapHelixFunctor::operator()(const Vector &Vect) const {
  return Vect.Z() - m_z;
}

ZPlaneHelixFunctor::ZPlaneHelixFunctor(double z,
				       const Vector &CellPosition,
				       bool Barrel):
  m_z(z),
  m_UnitZ(GetUnitZ(Barrel, CellPosition)),
  m_CellPositionZ(m_UnitZ.Dot(CellPosition)) {
}

double ZPlaneHelixFunctor::operator()(const Vector &Vect) const {
  return Vect.Dot(m_UnitZ) - m_CellPositionZ - m_z;
}

Vector ZPlaneHelixFunctor::GetUnitZ(bool Barrel, const Vector &Position) {
  if(Barrel) {
    const Vector PositionProjection(Position.X(), Position.Y(), 0.0);
    return PositionProjection.Unit();
  } else {
    return Vector(0.0, 0.0, 1.0);
  }
}

MirrorHelixFunctor::MirrorHelixFunctor(const Vector &MirrorPosition,
				       double Curvature):
  m_MirrorPosition(MirrorPosition),
  m_Curvature(Curvature) {
}

double MirrorHelixFunctor::operator()(const Vector &Vect) const {
  return (Vect - m_MirrorPosition).Mag2() - m_Curvature*m_Curvature;
}
