// Martin Duy Tat 28th October 2022

#include"TMath.h"
#include"Math/Functor.h"
#include"Math/RootFinder.h"
#include"Math/RotationZ.h"
#include"HelixPath.h"
#include"HelixFunctor.h"

using Functor1D = ROOT::Math::Functor1D;

HelixPath::HelixPath(const Vector &Momentum, int Q, double B):
  m_n(Momentum.Unit()),
  m_Radius(TMath::Sqrt(Momentum.Mag2())/(Q*B)),
  m_InvRadius(1.0/m_Radius),
  m_BFieldOn(B != 0.0),
  m_Origin_s(0.0) {
}

Vector HelixPath::GetPosition(double s) const {
  s -= m_Origin_s;
  const double Sin = TMath::Sin(s*m_InvRadius);
  const double Cos = TMath::Cos(s*m_InvRadius);
  const double x = m_BFieldOn ?
                   m_Radius*(m_n.X()*Sin + m_n.Y()*(1 - Cos)) :
                   m_n.X()*s;
  const double y = m_BFieldOn ?
                   m_Radius*(m_n.Y()*Sin - m_n.X()*(1 - Cos)) :
                   m_n.Y()*s;
  const double z = m_n.Z()*s;
  return m_Origin + Vector(x, y, z);
}

Vector HelixPath::GetDirection(double s) const {
  s -= m_Origin_s;
  const double Sin = TMath::Sin(s*m_InvRadius);
  const double Cos = TMath::Cos(s*m_InvRadius);
  const double x = m_BFieldOn ?
                   m_n.X()*Cos + m_n.Y()*Sin :
                   m_n.X();
  const double y = m_BFieldOn ?
                   m_n.Y()*Cos - m_n.X()*Sin :
                   m_n.Y();
  const double z = m_n.Z();
  return Vector(x, y, z);
}

double HelixPath::SolvePathLength(const HelixFunctor &Functor,
				  double Min, double Max) const {
  auto fcn = [&](double s) {
    return Functor(GetPosition(s));
  };
  ROOT::Math::RootFinder Solver;
  if(!Solver.Solve(fcn, Min, Max)) {
    return -999.0;
  }
  return Solver.Root();
}

void HelixPath::ResetOrigin(double NewOrigin_s) {
  m_Origin = GetPosition(NewOrigin_s);
  m_n = GetDirection(NewOrigin_s).Unit();
  m_Origin_s = NewOrigin_s;
}

void HelixPath::MapPhi(double DeltaPhi) {
  const RotationZ PhiRotation(DeltaPhi);
  m_n = PhiRotation(m_n);
  m_Origin = PhiRotation(m_Origin);
}

void HelixPath::ReflectZ() {
  m_n.SetZ(-m_n.Z());
  m_Origin.SetZ(-m_Origin.Z());
}

void HelixPath::ReflectY() {
  m_n.SetY(-m_n.Y());
  m_Origin.SetY(-m_Origin.Y());
}
