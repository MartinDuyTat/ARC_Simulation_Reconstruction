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
#include"Photon.h"

RadiatorCell::RadiatorCell(int CellNumber): m_ThetaLength(Settings::GetDouble("ARCGeometry/Length")/Settings::GetInt("ARCGeometry/ThetaCells")),
					    m_RadiatorThickness(Settings::GetDouble("RadiatorCell/RadiatorThickness")),
					    m_VesselThickness(Settings::GetDouble("RadiatorCell/VesselThickness")),
					    m_CoolingThickness(Settings::GetDouble("RadiatorCell/CoolingThickness")),
					    m_AerogelThickness(Settings::GetDouble("RadiatorCell/AerogelThickness")),
					    m_MirrorCurvature(Settings::GetDouble("RadiatorCell/MirrorCurvature")),
					    m_MirrorCentre(0.0, 0.0, GetMirrorCurvatureCentreZ()),
					    m_Position(GetCellPosition(CellNumber)),
                                            m_DeltaPhi(2.0*TMath::Pi()/Settings::GetInt("ARCGeometry/PhiCells")),
                                            m_CellNumber(CellNumber) {
  // TODO: Allow for off-axis mirror or mirror with different radius of curvature
  const int NumberThetaCells = Settings::GetInt("ARCGeometry/ThetaCells");
  if(m_CellNumber == -NumberThetaCells/2) {
    m_FirstMiddleLast = FirstMiddleLast::First;
  } else if(m_CellNumber == NumberThetaCells/2) {
    m_FirstMiddleLast = FirstMiddleLast::Last;
  } else {
    m_FirstMiddleLast = FirstMiddleLast::Middle;
  }
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

bool RadiatorCell::IsInsideThetaBoundary(const Vector &Position) const {
  return TMath::Abs(Position.X()) <= m_ThetaLength/2.0;
}

bool RadiatorCell::IsInsideThetaBoundary(const Photon &photon) const {
  return IsInsideThetaBoundary(photon.m_Position);
}

bool RadiatorCell::IsInsidePhiBoundary(const Photon &photon) const {
  const double PhiLength = m_DeltaPhi*Settings::GetDouble("ARCGeometry/Radius");
  return TMath::Abs(photon.m_Position.Y()) <= PhiLength/2.0 + photon.m_Position.Z()*TMath::Tan(m_DeltaPhi/2.0);
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>> RadiatorCell::DrawRadiatorGeometry() const {
  const auto MirrorCentreGlobal = m_Position + m_MirrorCentre;
  const double ArcAngle = TMath::ASin(0.5*m_ThetaLength/m_MirrorCurvature)*180.0/TMath::Pi();
  TArc MirrorArc(MirrorCentreGlobal.X(),
		 MirrorCentreGlobal.Z(),
		 m_MirrorCurvature,
		 90.0 - ArcAngle, 90.0 + ArcAngle);
  MirrorArc.SetLineColor(kBlack);
  MirrorArc.SetLineWidth(2);
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> Objects;
  Objects.push_back(std::make_pair(std::make_unique<TArc>(MirrorArc), "ONLY"));
  TLine LeftLine(-m_ThetaLength/2.0 + m_Position.X(),
		 Settings::GetDouble("ARCGeometry/Radius"),
		 -m_ThetaLength/2.0 + m_Position.X(),
		 Settings::GetDouble("ARCGeometry/Radius") + m_RadiatorThickness);
  LeftLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(LeftLine), ""));
  TLine RightLine(m_ThetaLength/2.0 + m_Position.X(),
		  Settings::GetDouble("ARCGeometry/Radius"),
		  m_ThetaLength/2.0 + m_Position.X(),
		  Settings::GetDouble("ARCGeometry/Radius") + m_RadiatorThickness);
  RightLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(RightLine), ""));
  TLine DetectorLine(-m_ThetaLength/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_VesselThickness + m_CoolingThickness,
		     m_ThetaLength/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_VesselThickness + m_CoolingThickness);
  DetectorLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(DetectorLine), ""));
  TLine AerogelLine(-m_ThetaLength/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_VesselThickness + m_CoolingThickness + m_AerogelThickness,
		     m_ThetaLength/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_VesselThickness + m_CoolingThickness + m_AerogelThickness);
  AerogelLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(AerogelLine), ""));
  return Objects;
}

Vector RadiatorCell::GetCellPosition(int CellNumber) const {
  const double ZPosition = Settings::GetDouble("ARCGeometry/Radius") + m_VesselThickness + m_CoolingThickness;
  if(CellNumber == 0) {
    return Vector(0.0, 0.0, ZPosition);
  } else {
    if(CellNumber > 0) {
      const double XPosition = m_ThetaLength*(CellNumber - 0.5);
      return Vector(XPosition, 0.0, ZPosition);
    } else {
      const double XPosition = m_ThetaLength*(CellNumber + 0.5);
      return Vector(XPosition, 0.0, ZPosition);
    }
  }
}

double RadiatorCell::GetThetaLength() const {
  return m_ThetaLength;
}

double RadiatorCell::GetCellNumber() const {
  return m_CellNumber;
}

RadiatorCell::FirstMiddleLast RadiatorCell::GetFirstMiddleLast() const {
  return m_FirstMiddleLast;
}
