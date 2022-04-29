// Martin Duy Tat 29th April 2022

#include"TrackingVolume.h"

TrackingVolume::TrackingVolume(double Radius, double FieldStrength): m_Radius(Radius),
								     m_FieldStrength(FieldStrength) {
}

double TrackingVolume::GetRadius() const {
  return m_Radius;
}

double TrackingVolume::GetFieldStrength() const {
  return m_FieldStrength;
}
