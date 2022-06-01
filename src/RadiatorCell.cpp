// Martin Duy Tat 29th April 2022

#include<vector>
#include<utility>
#include<string>
#include<memory>
#include"TArc.h"
#include"TLine.h"
#include"TMath.h"
#include"TObject.h"
#include"RadiatorCell.h"
#include"Settings.h"

RadiatorCell::RadiatorCell(const Vector &Position): m_RadiatorThickness(Settings::GetDouble("RadiatorCell/RadiatorThickness")),
						    m_VesselThickness(Settings::GetDouble("RadiatorCell/VesselThickness")),
						    m_CoolingThickness(Settings::GetDouble("RadiatorCell/CoolingThickness")),
						    m_AerogelThickness(Settings::GetDouble("RadiatorCell/AerogelThickness")),
						    m_MirrorCurvature(Settings::GetDouble("RadiatorCell/MirrorCurvature")),
						    m_MirrorCentre(0.0, 0.0, GetMirrorCurvatureCentreZ()),
						    m_Position(Position) {
  // TODO: Allow for off-axis mirror
}

double RadiatorCell::GetRadiatorThickness() const {
  return m_RadiatorThickness;
}

double RadiatorCell::GetVesselThickness() const {
  return m_VesselThickness;
}

double RadiatorCell::GetCoolingThickness() const {
  return m_CoolingThickness;
}

double RadiatorCell::GetAerogelThickness() const {
  return m_AerogelThickness;
}

const Vector& RadiatorCell::GetMirrorCentre() const {
  return m_MirrorCentre;
}

double RadiatorCell::GetMirrorCurvatureCentreZ() const {
  // TODO: Allow for off-axis mirror
  return m_RadiatorThickness - m_MirrorCurvature - 2*m_VesselThickness - m_CoolingThickness;
}

double RadiatorCell::GetMirrorCurvature() const {
  return m_MirrorCurvature;
}

const Vector& RadiatorCell::GetRadiatorPosition() const {
  return m_Position;
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>> RadiatorCell::DrawRadiatorGeometry() const {
  const auto MirrorCentreGlobal = m_Position + m_MirrorCentre;
  const double ARCLength = Settings::GetDouble("ARCGeometry/Length");
  const int ThetaCells = Settings::GetInt("ARCGeometry/ThetaCells");
  const double ThetaLength = ARCLength/static_cast<double>(ThetaCells);
  const double ArcAngle = TMath::ASin(0.5*ThetaLength/m_MirrorCurvature)*180.0/TMath::Pi();
  TArc MirrorArc(MirrorCentreGlobal.Y(),
		 MirrorCentreGlobal.Z(),
		 m_MirrorCurvature,
		 90.0 - ArcAngle, 90.0 + ArcAngle);
  MirrorArc.SetLineColor(kBlack);
  MirrorArc.SetLineWidth(2);
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> Objects;
  Objects.push_back(std::make_pair(std::make_unique<TArc>(MirrorArc), "ONLY"));
  TLine LeftLine(-ThetaLength/2.0,
		 Settings::GetDouble("ARCGeometry/Radius"),
		 -ThetaLength/2.0,
		 Settings::GetDouble("ARCGeometry/Radius") + m_RadiatorThickness);
  LeftLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(LeftLine), ""));
  TLine RightLine(ThetaLength/2.0,
		  Settings::GetDouble("ARCGeometry/Radius"),
		  ThetaLength/2.0,
		  Settings::GetDouble("ARCGeometry/Radius") + m_RadiatorThickness);
  RightLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(RightLine), ""));
  return Objects;
}
