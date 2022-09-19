// Martin Duy Tat 6th September 2022

#include<vector>
#include<memory>
#include<utility>
#include<string>
#include"TArc.h"
#include"TLine.h"
#include"TMath.h"
#include"TObject.h"
#include"BarrelRadiatorCell.h"
#include"RadiatorCell.h"
#include"Settings.h"

BarrelRadiatorCell::BarrelRadiatorCell(std::size_t CellColumnNumber,
				       std::size_t CellRowNumber,
				       double HexagonSize):
  RadiatorCell(CellColumnNumber, CellRowNumber, HexagonSize),
  m_Position(GetCellPosition(CellColumnNumber, CellRowNumber)) {
}

const Vector& BarrelRadiatorCell::GetRadiatorPosition() const {
  return m_Position;
}

bool BarrelRadiatorCell::IsInsideCell(const Vector &Position) const {
  // Get x and y coordinates after mapping everything to first quadrant
  const double x = TMath::Abs(Position.X());
  // Need a stretching factor in y direction because of the curvature
  /*const double Radius = m_Position.Z() - m_CoolingThickness;
  const double TanTheta = Position.Y()/(Radius + m_CoolingThickness);
  const double SecTheta = TMath::Sqrt(1 + TanTheta*TanTheta);
  const double Stretch = SecTheta*(1 + m_CoolingThickness/Radius);
  const double y = TMath::Abs(Position.Y())/Stretch;*/
  const double Radius = m_Position.Z() + Position.Z();
  const double Theta = TMath::Abs(Position.Y())/Radius;
  const double y = Radius*TMath::Sin(Theta);
  // First part is checking the sloped part, the other is the vertical part
  return x < std::min(m_HexagonSize - y*TMath::Sqrt(3.0), m_HexagonSize*0.5);
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
BarrelRadiatorCell::DrawRadiatorGeometry() const {
  // First find the mirror intersections with the walls in local coordinates
  // by solving a quadratic
  const double b = m_MirrorCentre.Z();
  const double c_left = TMath::Power(m_HexagonSize/2.0, 2)
                      + m_MirrorCentre.Mag2()
                      - m_MirrorCurvature*m_MirrorCurvature
                      + m_MirrorCentre.X()*m_HexagonSize;
  const double c_right = TMath::Power(m_HexagonSize/2.0, 2)
                       + m_MirrorCentre.Mag2()
                       - m_MirrorCurvature*m_MirrorCurvature
                       - m_MirrorCentre.X()*m_HexagonSize;
  const double s_left = TMath::Sqrt(b*b - c_left);
  const double s_right = TMath::Sqrt(b*b - c_right);
  const double SinLeft = s_left/m_MirrorCurvature;
  const double SinRight = s_right/m_MirrorCurvature;
  double LeftAngle = m_MirrorCentre.X() < -m_HexagonSize/2.0 ?
		     TMath::ASin(SinLeft)*180.0/TMath::Pi() :
		     180.0 - TMath::ASin(SinLeft)*180.0/TMath::Pi();
  double RightAngle = m_MirrorCentre.X() > m_HexagonSize/2.0 ?
		      180.0 - TMath::ASin(SinRight)*180.0/TMath::Pi() :
		      TMath::ASin(SinRight)*180.0/TMath::Pi();
  // Check if mirror intersects with the top wall
  const double bb = m_MirrorCentre.X();
  const double h = m_RadiatorThickness - 2.0*m_VesselThickness - m_CoolingThickness;
  const double cc = h*h + m_MirrorCentre.Mag2() - 2.0*m_MirrorCentre.Z()*h
                  - m_MirrorCurvature*m_MirrorCurvature;
  const double s1 = bb + TMath::Sqrt(bb*bb - cc);
  const double s2 = bb - TMath::Sqrt(bb*bb - cc);
  if(bb*bb > cc &&
     (TMath::Abs(s1) < m_HexagonSize/2.0 || TMath::Abs(s2) < m_HexagonSize/2.0)) {
    const double s = TMath::Abs(s1) < m_HexagonSize/2.0 ? s1 : s2;
    const double MidAngle = TMath::ATan2(h - m_MirrorCentre.Z(),
					 s - m_MirrorCentre.X())*180.0/TMath::Pi();
    if(m_MirrorCentre.X() < 0.0) {
      LeftAngle = MidAngle;
    } else {
      RightAngle = MidAngle;
    }
  }
  // Finally draw everything
  const auto MirrorCentreGlobal = m_Position + m_MirrorCentre;
  TArc MirrorArc(MirrorCentreGlobal.X(),
		 MirrorCentreGlobal.Z(),
		 m_MirrorCurvature,
		 RightAngle,
		 LeftAngle);
  MirrorArc.SetLineColor(kBlack);
  MirrorArc.SetLineWidth(2);
  MirrorArc.SetFillColorAlpha(kWhite, 0.0);
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> Objects;
  Objects.push_back(std::make_pair(std::make_unique<TArc>(MirrorArc), "ONLY"));
  // Draw the walls
  TLine LeftLine(-m_HexagonSize/2.0 + m_Position.X(),
		 Settings::GetDouble("ARCGeometry/Radius") - m_VesselThickness,
		 -m_HexagonSize/2.0 + m_Position.X(),
		 Settings::GetDouble("ARCGeometry/Radius")
		 + m_RadiatorThickness - m_VesselThickness);
  LeftLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(LeftLine), ""));
  TLine RightLine(m_HexagonSize/2.0 + m_Position.X(),
		  Settings::GetDouble("ARCGeometry/Radius") - m_VesselThickness,
		  m_HexagonSize/2.0 + m_Position.X(),
		  Settings::GetDouble("ARCGeometry/Radius")
		  + m_RadiatorThickness - m_VesselThickness);
  RightLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(RightLine), ""));
  // Draw the detector plane
  TLine DetectorLine(-m_HexagonSize/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_CoolingThickness,
		     m_HexagonSize/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_CoolingThickness);
  DetectorLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(DetectorLine), ""));
  // Draw the aerogel plane
  TLine AerogelLine(-m_HexagonSize/2.0 + m_Position.X(),
		    Settings::GetDouble("ARCGeometry/Radius")
		    + m_CoolingThickness + m_AerogelThickness,
		    m_HexagonSize/2.0 + m_Position.X(),
		    Settings::GetDouble("ARCGeometry/Radius")
		    + m_CoolingThickness + m_AerogelThickness);
  AerogelLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(AerogelLine), ""));
  // Draw the bottom cooling plane
  TLine CoolingLine(-m_HexagonSize/2.0 + m_Position.X(),
		    Settings::GetDouble("ARCGeometry/Radius"),
		    m_HexagonSize/2.0 + m_Position.X(),
		    Settings::GetDouble("ARCGeometry/Radius"));
  CoolingLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(CoolingLine), ""));
  // Draw SiPM
  Objects.push_back(std::make_pair(m_Detector.DrawSiPM(GetRadiatorPosition()), ""));
  return Objects;
}

Vector BarrelRadiatorCell::GetCellPosition(std::size_t CellColumnNumber,
					   std::size_t CellRowNumber) const {
  if(CellColumnNumber > 9) {
    throw std::invalid_argument("Invalid cell column number: "
				+ std::to_string(CellColumnNumber));
  }
  const double ZPosition = Settings::GetDouble("ARCGeometry/Radius")
                         + m_CoolingThickness;
  if(CellRowNumber == 0 && CellColumnNumber == 0) {
    return Vector(0.0, 0.0, ZPosition);
  }
  if(CellRowNumber == 1) {
    const double XPosition = m_HexagonSize*static_cast<double>(CellColumnNumber);
    return Vector(XPosition, 0.0, ZPosition);
  } else if(CellRowNumber == 2) {
    const double XPosition = m_HexagonSize*(static_cast<double>(CellColumnNumber) - 0.5);
    return Vector(XPosition, 0.0, ZPosition);
  } else {
    throw std::invalid_argument("Invalid cell row number: "
				+ std::to_string(CellRowNumber));
  }
}
