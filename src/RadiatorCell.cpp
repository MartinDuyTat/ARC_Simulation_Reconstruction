// Martin Duy Tat 29th April 2022

#include"RadiatorCell.h"

RadiatorCell::RadiatorCell(const Vector &Position): m_RadiatorThickness(0.20),
						    m_VesselThickness(0.02),
						    m_CoolingThickness(0.02),
						    m_AerogelThickness(0.02),
						    m_MirrorCurvature(0.30),
						    m_MirrorCentre(0.0, 0.0, GetMirrorCurvatureCentreZ()),
						    m_Position(Position) {
  // TODO: Allow for off-axis mirror
}

double RadiatorCell::GetRadiatorThickness() const {
  return m_RadiatorThickness;
}

double VesselCell::GetVesselThickness() const {
  return m_VesselThickness;
}

double CoolingCell::GetCoolingThickness() const {
  return m_CoolingThickness;
}

double AerogelCell::GetAerogelThickness() const {
  return m_AerogelThickness;
}

double RadiatorCell::GetMirrorCurvatureCentreZ() const {
  return m_RadiatorThickness - m_MirrorCurvature - 2*VesselThickness - m_CoolingThickness;
}

const Vector& RadiatorCell::GetRadiatorPosition() const {
  return m_Position;
}
