// Martin Duy Tat 1st September 2022

#include<array>
#include<utility>
#include"TMath.h"
#include"TArc.h"
#include"EndCapRadiatorCell.h"
#include"RadiatorCell.h"
#include"Settings.h"
#include"BarrelRadiatorCell.h"
#include"Utilities.h"

EndCapRadiatorCell::EndCapRadiatorCell(int CellColumnNumber,
				       int CellRowNumber,
				       double HexagonSize):
  RadiatorCell::RadiatorCell(CellColumnNumber, CellRowNumber, HexagonSize, "EndCap"),
  m_HexagonDistY(m_HexagonSize*TMath::Sqrt(3)/2.0),
  m_Position(GetCellPosition(CellColumnNumber, CellRowNumber)) {
  const std::string RadiatorName = "EndCapRadiator_c"
                                 + std::to_string(m_CellNumber.first)
                                 + "_r"
                                 + std::to_string(m_CellNumber.second)
                                 + "_";
  const std::string XPositionName = "RadiatorCell/" + RadiatorName + "XPosition";
  if(Settings::Exists(XPositionName)) {
    SetMirrorXPosition(Settings::GetDouble(XPositionName));
  }
}

bool EndCapRadiatorCell::IsInsideCell(const Vector &Position) const {
  // Get x and y coordinates after mapping everything to first quadrant
  const double x = TMath::Abs(Position.X());
  const double y = TMath::Abs(Position.Y());
  // First part is checking the sloped part, the other is the vertical part
  return x < std::min(m_HexagonSize - y*TMath::Sqrt(3.0), m_HexagonSize*0.5);
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
EndCapRadiatorCell::DrawRadiatorGeometry() const {
  // First find the mirror intersections with the walls in local coordinates
  // by solving a quadratic
  const double b = m_MirrorCentre.Z();
  const double c_left = TMath::Power(m_HexagonSize/2.0, 2)
                      + m_MirrorCentre.Mag2()
                      - m_MirrorCurvature*m_MirrorCurvature
                      - m_MirrorCentre.X()*m_HexagonSize;
  const double c_right = TMath::Power(m_HexagonSize/2.0, 2)
                       + m_MirrorCentre.Mag2()
                       - m_MirrorCurvature*m_MirrorCurvature
                       + m_MirrorCentre.X()*m_HexagonSize;
  const double s_left = TMath::Sqrt(b*b - c_left);
  const double s_right = TMath::Sqrt(b*b - c_right);
  const double SinLeft = s_left/m_MirrorCurvature;
  const double SinRight = s_right/m_MirrorCurvature;
  double LeftAngle = -m_MirrorCentre.X() < -m_HexagonSize/2.0 ?
		     TMath::ASin(SinLeft)*180.0/TMath::Pi() :
		     180.0 - TMath::ASin(SinLeft)*180.0/TMath::Pi();
  double RightAngle = -m_MirrorCentre.X() > m_HexagonSize/2.0 ?
		      180.0 - TMath::ASin(SinRight)*180.0/TMath::Pi() :
		      TMath::ASin(SinRight)*180.0/TMath::Pi();
  // Check if mirror intersects with the top wall
  const double bb = -m_MirrorCentre.X();
  const double h = m_RadiatorThickness - 2.0*m_VesselThickness - m_CoolingThickness;
  const double cc = h*h + m_MirrorCentre.Mag2() - 2.0*m_MirrorCentre.Z()*h
                  - m_MirrorCurvature*m_MirrorCurvature;
  const double s1 = bb + TMath::Sqrt(bb*bb - cc);
  const double s2 = bb - TMath::Sqrt(bb*bb - cc);
  if(bb*bb > cc &&
     (TMath::Abs(s1) < m_HexagonSize/2.0 || TMath::Abs(s2) < m_HexagonSize/2.0)) {
    const double s = TMath::Abs(s1) < m_HexagonSize/2.0 ? s1 : s2;
    const double MidAngle = TMath::ATan2(h - m_MirrorCentre.Z(),
					 s + m_MirrorCentre.X())*180.0/TMath::Pi();
    if(m_MirrorCentre.X() > 0.0) {
      LeftAngle = MidAngle;
    } else {
      RightAngle = MidAngle;
    }
  }
  // Finally draw everything
  const auto MirrorCentreGlobal = Utilities::SwapXZForEndCap(this,
							     m_Position
							     + m_MirrorCentre);
  const auto SwappedPosition = Utilities::SwapXZForEndCap(this, m_Position);
  TArc MirrorArc(MirrorCentreGlobal.X(),
		 MirrorCentreGlobal.Z(),
		 m_MirrorCurvature,
		 RightAngle - 90.0,
		 LeftAngle - 90.0);
  MirrorArc.SetLineColor(kBlack);
  MirrorArc.SetLineWidth(2);
  MirrorArc.SetFillColorAlpha(kWhite, 0.0);
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> Objects;
  Objects.push_back(std::make_pair(std::make_unique<TArc>(MirrorArc), "ONLY"));
  // Draw the walls
  TLine LeftLine(Settings::GetDouble("ARCGeometry/BarrelZ") - m_VesselThickness,
		 -m_HexagonSize/2.0 + SwappedPosition.Z(),
		 Settings::GetDouble("ARCGeometry/BarrelZ")
		 + m_RadiatorThickness - m_VesselThickness,
		 -m_HexagonSize/2.0 + SwappedPosition.Z());
  LeftLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(LeftLine), ""));
  TLine RightLine(Settings::GetDouble("ARCGeometry/BarrelZ") - m_VesselThickness,
		  m_HexagonSize/2.0 + SwappedPosition.Z(),
		  Settings::GetDouble("ARCGeometry/BarrelZ")
		  + m_RadiatorThickness - m_VesselThickness,
		  m_HexagonSize/2.0 + SwappedPosition.Z());
  RightLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(RightLine), ""));
  // Draw the detector plane
  TLine DetectorLine(Settings::GetDouble("ARCGeometry/BarrelZ") + m_CoolingThickness,
		     -m_HexagonSize/2.0 + SwappedPosition.Z(),
		     Settings::GetDouble("ARCGeometry/BarrelZ") + m_CoolingThickness,
		     m_HexagonSize/2.0 + SwappedPosition.Z());
  DetectorLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(DetectorLine), ""));
  // Draw the aerogel plane
  TLine AerogelLine(Settings::GetDouble("ARCGeometry/BarrelZ")
		    + m_CoolingThickness + m_AerogelThickness,
		    -m_HexagonSize/2.0 + SwappedPosition.Z(),
		    Settings::GetDouble("ARCGeometry/BarrelZ")
		    + m_CoolingThickness + m_AerogelThickness,
		    m_HexagonSize/2.0 + SwappedPosition.Z());
  AerogelLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(AerogelLine), ""));
  // Draw the bottom cooling plane
  TLine CoolingLine(Settings::GetDouble("ARCGeometry/BarrelZ"),
		    -m_HexagonSize/2.0 + SwappedPosition.Z(),
		    Settings::GetDouble("ARCGeometry/BarrelZ"),
		    m_HexagonSize/2.0 + SwappedPosition.Z());
  CoolingLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(CoolingLine), ""));
  // Draw SiPM
  auto DrawnDetector = m_Detector.DrawSiPM(GetRadiatorPosition());
  Objects.push_back(std::make_pair(std::move(DrawnDetector), ""));
  return Objects;
}

const Vector& EndCapRadiatorCell::GetRadiatorPosition() const {
  return m_Position;
}

Vector EndCapRadiatorCell::GetCellPosition(int CellColumnNumber,
					   int CellRowNumber) const {
  if(std::find(m_ValidCells.begin(),
	       m_ValidCells.end(),
	       std::make_pair(CellColumnNumber, CellRowNumber)) ==
     m_ValidCells.end()) {
    throw std::invalid_argument("Invalid cell column number: "
				+ std::to_string(CellColumnNumber));
  }
  const double ZPosition = Settings::GetDouble("ARCGeometry/BarrelZ")
                         + m_CoolingThickness;
  const double YPosition = m_HexagonDistY*(CellRowNumber - 1);
  if(CellRowNumber%2 == 1) {
    const double XPosition = m_HexagonSize*CellColumnNumber;
    return Vector(XPosition, YPosition, ZPosition);
  } else {
    const double XPosition = m_HexagonSize*(CellColumnNumber - 0.5);
    return Vector(XPosition, YPosition, ZPosition);
  }
}

void EndCapRadiatorCell::SetMirrorXPosition(double x) {
  m_MirrorCentre.SetX(m_DefaultMirrorCentre.X() - x);
}
