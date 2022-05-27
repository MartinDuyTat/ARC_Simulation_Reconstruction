// Martin Duy Tat 29th April 2022

#include"TObject.h"
#include"TLine.h"
#include"TrackingVolume.h"
#include"Settings.h"

TrackingVolume::TrackingVolume(): m_Radius(Settings::GetDouble("ARCGeometry/Radius")),
				  m_Length(Settings::GetDouble("ARCGeometry/Length")),
				  m_ThetaCells(Settings::GetInt("ARCGeometry/ThetaCells")),
				  m_PhiCells(Settings::GetInt("ARCGeometry/PhiCells")),
				  m_FieldStrength(Settings::GetDouble("ARCGeometry/FieldStrength")) {
}

double TrackingVolume::GetRadius() const {
  return m_Radius;
}

double TrackingVolume::GetFieldStrength() const {
  return m_FieldStrength;
}

std::vector<std::unique_ptr<TObject>> TrackingVolume::DrawARCGeometry() const {
  std::vector<std::unique_ptr<TObject>> Lines;
  const double ThetaLength = m_Length/m_ThetaCells;
  TLine Line1(0.0, 0.0, ThetaLength/2.0, m_Radius);
  Line1.SetLineStyle(kDashed);
  Line1.SetLineWidth(1);
  Line1.SetLineColor(kBlack);
  Lines.push_back(std::make_unique<TLine>(Line1));
  TLine Line2(0.0, 0.0, -ThetaLength/2.0, m_Radius);
  Line2.SetLineStyle(kDashed);
  Line2.SetLineWidth(1);
  Line2.SetLineColor(kBlack);
  Lines.push_back(std::make_unique<TLine>(Line2));
  return Lines;
}
