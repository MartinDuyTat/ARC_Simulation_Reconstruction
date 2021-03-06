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
					    m_Position(GetCellPosition(CellNumber)),
					    m_Detector(DetermineSiPMPositionX(), 0.0),
					    m_MirrorCurvature(DetermineMirrorCurvature()),
					    m_MirrorCentre(0.0, 0.0, GetMirrorCurvatureCentreZ()),
                                            m_DeltaPhi(2.0*TMath::Pi()/Settings::GetInt("ARCGeometry/PhiCells")),
                                            m_CellNumber(CellNumber) {
  // TODO: Allow for off-axis mirror or mirror with different radius of curvature
  // Check if this cell is at the edge
  if(Settings::GetBool("General/FullArray")) {
    const int NumberThetaCells = Settings::GetInt("ARCGeometry/ThetaCells");
    if(m_CellNumber == -NumberThetaCells/2) {
      m_FirstMiddleLast = FirstMiddleLast::First;
    } else if(m_CellNumber == NumberThetaCells/2) {
      m_FirstMiddleLast = FirstMiddleLast::Last;
    } else {
      m_FirstMiddleLast = FirstMiddleLast::Middle;
    }
  } else {
    m_FirstMiddleLast = FirstMiddleLast::Single;
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
  // First find the mirror intersections with the walls in local coordinates by solving a quadratic
  const double b = m_MirrorCentre.Z();
  const double c_left = TMath::Power(m_ThetaLength/2.0, 2)
                      + m_MirrorCentre.Mag2()
                      - m_MirrorCurvature*m_MirrorCurvature
                      - m_MirrorCentre.X()*m_ThetaLength;
  const double c_right = TMath::Power(m_ThetaLength/2.0, 2)
                       + m_MirrorCentre.Mag2()
                       - m_MirrorCurvature*m_MirrorCurvature
                       + m_MirrorCentre.X()*m_ThetaLength;
  const double s_left = TMath::Sqrt(b*b - c_left);
  const double s_right = TMath::Sqrt(b*b - c_right);
  const auto MirrorCentreGlobal = m_Position + m_MirrorCentre;
  //const double ArcAngle = TMath::ASin(0.5*m_ThetaLength/m_MirrorCurvature)*180.0/TMath::Pi();
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
  // Draw the detector plane
  TLine DetectorLine(-m_ThetaLength/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_VesselThickness + m_CoolingThickness,
		     m_ThetaLength/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_VesselThickness + m_CoolingThickness);
  DetectorLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(DetectorLine), ""));
  // Draw the aerogel plane
  TLine AerogelLine(-m_ThetaLength/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_VesselThickness + m_CoolingThickness + m_AerogelThickness,
		     m_ThetaLength/2.0 + m_Position.X(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_VesselThickness + m_CoolingThickness + m_AerogelThickness);
  AerogelLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(AerogelLine), ""));
  // Draw SiPM
  Objects.push_back(std::make_pair(m_Detector.DrawSiPM(GetRadiatorPosition()), ""));
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

double RadiatorCell::DetermineMirrorCurvature() const {
  const double NominalCurvature = Settings::GetDouble("RadiatorCell/MirrorCurvature");
  if(Settings::GetBool("General/VariableMirrorCurvature")) {
    const double Radius = Settings::GetDouble("ARCGeometry/Radius");
    const double Theta = TMath::ATan2(m_Position.X(), Radius + m_RadiatorThickness - m_VesselThickness);
    return NominalCurvature/TMath::Abs(TMath::Cos(Theta));
  } else {
    return NominalCurvature;
  }
}

double RadiatorCell::DetermineSiPMPositionX() const {
  if(Settings::GetBool("General/VariableMirrorCurvature")) {
    const double TanTheta = GetRadiatorPosition().X()/GetRadiatorPosition().Z();
    const double SiPM_xPosition = 2*(m_RadiatorThickness - 2*m_VesselThickness - m_CoolingThickness)*TanTheta;
    return SiPM_xPosition;
  } else {
    return 0.0;
  }
}

SiPM& RadiatorCell::GetDetector() {
  return m_Detector;
}
