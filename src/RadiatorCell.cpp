// Martin Duy Tat 29th April 2022

#include<vector>
#include<utility>
#include<string>
#include<memory>
#include<stdexcept>
#include"TArc.h"
#include"TLine.h"
#include"TMath.h"
#include"TObject.h"
#include"RadiatorCell.h"
#include"Settings.h"
#include"Photon.h"

RadiatorCell::RadiatorCell(int CellColumnNumber,
			   int CellRowNumber,
			   double HexagonSize):
			   m_RadiatorThickness(Settings::GetDouble("RadiatorCell/RadiatorThickness")),
			   m_VesselThickness(Settings::GetDouble("RadiatorCell/VesselThickness")),
			   m_CoolingThickness(Settings::GetDouble("RadiatorCell/CoolingThickness")),
			   m_AerogelThickness(Settings::GetDouble("RadiatorCell/AerogelThickness")),
			   m_HexagonSize(HexagonSize),
			   m_CellNumber(std::make_pair(CellColumnNumber, CellRowNumber)),
			   m_Position(GetCellPosition(CellColumnNumber, CellRowNumber)),
			   m_Detector(DetermineSiPMPositionX(), 0.0),
			   m_MirrorCurvature(DetermineMirrorCurvature()),
			   m_MirrorCentre(0.0, 0.0, GetMirrorCurvatureCentreZ()) {
  // TODO: Allow for off-axis mirror or mirror with different radius of curvature
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

bool RadiatorCell::IsInsideCell(const Vector &Position) const {
  // Get x and y coordinates after mapping everything to first quadrant
  const double x = TMath::Abs(Position.X());
  // Need a stretching factor in y direction because of the curvature
  const double Radius = m_Position.Z() - m_CoolingThickness;
  const double TanTheta = Position.Y()/(Radius + m_CoolingThickness);
  const double SecTheta = TMath::Sqrt(1 + TanTheta*TanTheta);
  const double Stretch = SecTheta*(1 + m_CoolingThickness/Radius);
  const double y = TMath::Abs(Position.Y())/Stretch;
  // First part is checking the sloped part, the other is the vertical part
  return x < std::min(m_HexagonSize - y*TMath::Sqrt(3.0), m_HexagonSize*0.5);
}

bool RadiatorCell::IsInsideCell(const Photon &photon) const {
  return IsInsideCell(photon.m_Position);
}


std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
RadiatorCell::DrawRadiatorGeometry() const {
  // First find the mirror intersections with the walls in local coordinates by solving a quadratic
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
  const auto MirrorCentreGlobal = m_Position + m_MirrorCentre;
  //const double ArcAngle = TMath::ASin(0.5*m_HexagonSize/m_MirrorCurvature)*180.0/TMath::Pi();
  TArc MirrorArc(MirrorCentreGlobal.X(),
		 MirrorCentreGlobal.Z(),
		 m_MirrorCurvature,
		 TMath::ASin(s_right/m_MirrorCurvature)*180.0/TMath::Pi(),
		 180.0 - TMath::ASin(s_left/m_MirrorCurvature)*180.0/TMath::Pi());
  MirrorArc.SetLineColor(kBlack);
  MirrorArc.SetLineWidth(2);
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

Vector RadiatorCell::GetCellPosition(int CellColumnNumber, int CellRowNumber) const {
  if(CellColumnNumber > 9 || CellColumnNumber < 0) {
    throw std::invalid_argument("Invalid cell column number: " + std::to_string(CellColumnNumber));
  }
  const double ZPosition = Settings::GetDouble("ARCGeometry/Radius") + m_CoolingThickness;
  if(CellRowNumber == 0 && CellColumnNumber == 0) {
    return Vector(0.0, 0.0, ZPosition);
  }
  if(CellRowNumber == 1) {
    const double XPosition = m_HexagonSize*CellColumnNumber;
    return Vector(XPosition, 0.0, ZPosition);
  } else if(CellRowNumber == 2) {
    const double XPosition = m_HexagonSize*(CellColumnNumber - 0.5);
    return Vector(XPosition, 0.0, ZPosition);
  } else {
    throw std::invalid_argument("Invalid cell row number: " + std::to_string(CellRowNumber));
  }
}

double RadiatorCell::GetHexagonSize() const {
  return m_HexagonSize;
}

std::pair<int, int> RadiatorCell::GetCellNumber() const {
  return m_CellNumber;
}

double RadiatorCell::DetermineMirrorCurvature() const {
  const double NominalCurvature = Settings::GetDouble("RadiatorCell/MirrorCurvature");
  return NominalCurvature;
  /*if(Settings::GetBool("General/VariableMirrorCurvature")) {
    const double Radius = Settings::GetDouble("ARCGeometry/Radius");
    const double Theta = TMath::ATan2(m_Position.X(), Radius + m_RadiatorThickness - m_VesselThickness);
    return NominalCurvature/TMath::Abs(TMath::Cos(Theta));
  } else {
    return NominalCurvature;
  }*/
}

double RadiatorCell::DetermineSiPMPositionX() const {
  /*if(Settings::GetBool("General/VariableMirrorCurvature")) {
    const double TanTheta = GetRadiatorPosition().X()/GetRadiatorPosition().Z();
    const double SiPM_xPosition = 2*(m_RadiatorThickness - 2*m_VesselThickness - m_CoolingThickness)*TanTheta;
    return SiPM_xPosition;
  } else {
    return 0.0;
  }*/
  return 0.0;
}

SiPM& RadiatorCell::GetDetector() {
  return m_Detector;
}
