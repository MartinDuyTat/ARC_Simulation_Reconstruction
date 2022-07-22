// Martin Duy Tat 29th April 2022

#include<vector>
#include<string>
#include<memory>
#include<utility>
#include"TObject.h"
#include"TLine.h"
#include"TrackingVolume.h"
#include"Settings.h"

TrackingVolume::TrackingVolume():
  m_Radius(Settings::GetDouble("ARCGeometry/Radius")),
  m_Length(Settings::GetDouble("ARCGeometry/Length")),
  m_CellsPerRow(Settings::GetInt("ARCGeometry/CellsPerRow")),
  m_FieldStrength(Settings::GetDouble("ARCGeometry/FieldStrength")) {
}

double TrackingVolume::GetRadius() const {
  return m_Radius;
}

double TrackingVolume::GetFieldStrength() const {
  return m_FieldStrength;
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
TrackingVolume::DrawARCGeometry() const {
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> Lines;
  const double ThetaLength = m_Length/(2*m_CellsPerRow - 1);
  const double VesselThickness = Settings::GetDouble("RadiatorCell/VesselThickness");
  if(Settings::GetBool("General/FullArray")) {
    for(int i = 0; i <= m_CellsPerRow; i++) {
      TLine Line1(0.0, 0.0, ThetaLength*(i - 0.5), m_Radius - VesselThickness);
      Line1.SetLineStyle(kDotted);
      Line1.SetLineWidth(1);
      Line1.SetLineColor(kBlack);
      Lines.push_back(std::make_pair(std::make_unique<TLine>(Line1), ""));
    }
  } else {
    TLine Line1(0.0, 0.0, ThetaLength/2.0, m_Radius - VesselThickness);
    Line1.SetLineStyle(kDashed);
    Line1.SetLineWidth(1);
    Line1.SetLineColor(kBlack);
    Lines.push_back(std::make_pair(std::make_unique<TLine>(Line1), ""));
    TLine Line2(0.0, 0.0, -ThetaLength/2.0, m_Radius - VesselThickness);
    Line2.SetLineStyle(kDashed);
    Line2.SetLineWidth(1);
    Line2.SetLineColor(kBlack);
    Lines.push_back(std::make_pair(std::make_unique<TLine>(Line2), ""));
  }
  return Lines;
}
