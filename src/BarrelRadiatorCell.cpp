// Martin Duy Tat 6th September 2022

#include<vector>
#include<memory>
#include<utility>
#include<string>
#include"TArc.h"
#include"TLine.h"
#include"TMath.h"
#include"Math/RotationY.h"
#include"Math/RotationZ.h"
#include"TObject.h"
#include"BarrelRadiatorCell.h"
#include"RadiatorCell.h"
#include"Settings.h"

using RotationY = ROOT::Math::RotationY;
using RotationZ = ROOT::Math::RotationZ;

BarrelRadiatorCell::BarrelRadiatorCell(std::size_t CellColumnNumber,
				       std::size_t CellRowNumber,
				       double HexagonSize):
  RadiatorCell(CellColumnNumber,
	       CellRowNumber,
	       HexagonSize,
	       GetCellPosition(CellColumnNumber, CellRowNumber, HexagonSize),
	       GetCellOrientation(HexagonSize, CellRowNumber)),
  m_BarrelRadius(Settings::GetDouble("ARCGeometry/Radius")
	       + Settings::GetDouble("RadiatorCell/CoolingThickness")) {
}

bool BarrelRadiatorCell::IsInsideCell(const ARCVector &Position) const {
  // Use global position to determine the radius
  const auto GlobalPosition = Position.GlobalVector();
  const double Radius = TMath::Sqrt(GlobalPosition.X()*GlobalPosition.X()
                                  + GlobalPosition.Y()*GlobalPosition.Y());
  // Get x and y coordinates after mapping everything to first quadrant
  const double x = TMath::Abs(Position.LocalVector().X());
  // Need to project y coordinate down to the barrel arc length
  const double ProjectedYAngle = TMath::ASin(Position.LocalVector().Y()/Radius);
  const double y = TMath::Abs(m_BarrelRadius*ProjectedYAngle);
  // First part is checking the sloped part, the other is the vertical part
  return x < std::min(m_HexagonSize - y*TMath::Sqrt(3.0), m_HexagonSize*0.5);
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
BarrelRadiatorCell::DrawRadiatorGeometry() const {
  // First find the mirror intersections with the walls in local coordinates
  // by solving a quadratic
  const auto MirrorCentre = GetMirrorCentre();
  const double b = MirrorCentre.Z();
  const double c_left = TMath::Power(m_HexagonSize/2.0, 2)
                      + MirrorCentre.Mag2()
                      - m_MirrorCurvature*m_MirrorCurvature
                      + MirrorCentre.X()*m_HexagonSize;
  const double c_right = TMath::Power(m_HexagonSize/2.0, 2)
                       + MirrorCentre.Mag2()
                       - m_MirrorCurvature*m_MirrorCurvature
                       - MirrorCentre.X()*m_HexagonSize;
  const double s_left = TMath::Sqrt(b*b - c_left);
  const double s_right = TMath::Sqrt(b*b - c_right);
  const double SinLeft = s_left/m_MirrorCurvature;
  const double SinRight = s_right/m_MirrorCurvature;
  double LeftAngle = MirrorCentre.X() < -m_HexagonSize/2.0 ?
		     TMath::ASin(SinLeft)*180.0/TMath::Pi() :
		     180.0 - TMath::ASin(SinLeft)*180.0/TMath::Pi();
  double RightAngle = MirrorCentre.X() > m_HexagonSize/2.0 ?
		      180.0 - TMath::ASin(SinRight)*180.0/TMath::Pi() :
		      TMath::ASin(SinRight)*180.0/TMath::Pi();
  // Check if mirror intersects with the top wall
  const double bb = MirrorCentre.X();
  const double h = m_RadiatorThickness - 2.0*m_VesselThickness - m_CoolingThickness;
  const double cc = h*h + MirrorCentre.Mag2() - 2.0*MirrorCentre.Z()*h
                  - m_MirrorCurvature*m_MirrorCurvature;
  const double s1 = bb + TMath::Sqrt(bb*bb - cc);
  const double s2 = bb - TMath::Sqrt(bb*bb - cc);
  if(bb*bb > cc &&
     (TMath::Abs(s1) < m_HexagonSize/2.0 || TMath::Abs(s2) < m_HexagonSize/2.0)) {
    const double s = TMath::Abs(s1) < m_HexagonSize/2.0 ? s1 : s2;
    const double MidAngle = TMath::ATan2(h - MirrorCentre.Z(),
					 s - MirrorCentre.X())*180.0/TMath::Pi();
    if(MirrorCentre.X() < 0.0) {
      LeftAngle = MidAngle;
    } else {
      RightAngle = MidAngle;
    }
  }
  // Finally draw everything
  const auto ReversePhi = ReversePhiRotation();
  const auto MirrorCentreGlobal = ReversePhi(m_MirrorCentre.GlobalVector());
  TArc MirrorArc(MirrorCentreGlobal.Z(),
		 MirrorCentreGlobal.X(),
		 m_MirrorCurvature,
		 RightAngle,
		 LeftAngle);
  MirrorArc.SetLineColor(kBlack);
  MirrorArc.SetLineWidth(2);
  MirrorArc.SetFillColorAlpha(kWhite, 0.0);
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> Objects;
  Objects.push_back(std::make_pair(std::make_unique<TArc>(MirrorArc), "ONLY"));
  // Draw the walls
  TLine LeftLine(-m_HexagonSize/2.0 + m_Position.Z(),
		 Settings::GetDouble("ARCGeometry/Radius") - m_VesselThickness,
		 -m_HexagonSize/2.0 + m_Position.Z(),
		 Settings::GetDouble("ARCGeometry/Radius")
		 + m_RadiatorThickness - m_VesselThickness);
  LeftLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(LeftLine), ""));
  TLine RightLine(m_HexagonSize/2.0 + m_Position.Z(),
		  Settings::GetDouble("ARCGeometry/Radius") - m_VesselThickness,
		  m_HexagonSize/2.0 + m_Position.Z(),
		  Settings::GetDouble("ARCGeometry/Radius")
		  + m_RadiatorThickness - m_VesselThickness);
  RightLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(RightLine), ""));
  // Draw the detector plane
  TLine DetectorLine(-m_HexagonSize/2.0 + m_Position.Z(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_CoolingThickness,
		     m_HexagonSize/2.0 + m_Position.Z(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_CoolingThickness);
  DetectorLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(DetectorLine), ""));
  // Draw the aerogel plane
  TLine AerogelLine(-m_HexagonSize/2.0 + m_Position.Z(),
		    Settings::GetDouble("ARCGeometry/Radius")
		    + m_CoolingThickness + m_AerogelThickness,
		    m_HexagonSize/2.0 + m_Position.Z(),
		    Settings::GetDouble("ARCGeometry/Radius")
		    + m_CoolingThickness + m_AerogelThickness);
  AerogelLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(AerogelLine), ""));
  // Draw the bottom cooling plane
  TLine CoolingLine(-m_HexagonSize/2.0 + m_Position.Z(),
		    Settings::GetDouble("ARCGeometry/Radius"),
		    m_HexagonSize/2.0 + m_Position.Z(),
		    Settings::GetDouble("ARCGeometry/Radius"));
  CoolingLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(CoolingLine), ""));
  // Draw SiPM
  const auto RadiatorPosition = ReversePhi(GetRadiatorPosition());
  Objects.push_back(std::make_pair(m_Detector.DrawSiPM(RadiatorPosition), ""));
  return Objects;
}

RotationZ BarrelRadiatorCell::ReversePhiRotation() const {
  if(m_CellNumber.second == 1) {
    return RotationZ(0.0);
  } else {
    const double DeltaPhi = 0.5*TMath::Sqrt(3)*m_HexagonSize/m_BarrelRadius;
    return RotationZ(-DeltaPhi);
  }
}

bool BarrelRadiatorCell::IsDetectorInsideCell() const {
  const double Position = m_Detector.GetDetectorXPosition();
  const double Tilt = m_Detector.GetDetectorTilt();
  const double DetSize = m_Detector.GetDetectorSizeX();
  if(Position + DetSize*TMath::Cos(Tilt)/2.0 > m_HexagonSize/2.0) {
    return false;
  } else if(Position - DetSize*TMath::Cos(Tilt)/2.0 < -m_HexagonSize/2.0) {
    return false;
  } else {
    return true;
  }
}

Vector BarrelRadiatorCell::GetCellPosition(std::size_t CellColumnNumber,
					   std::size_t CellRowNumber,
					   double HexagonSize) {
  if(CellColumnNumber > 9) {
    throw std::invalid_argument("Invalid cell column number: "
				+ std::to_string(CellColumnNumber));
  }
  const double XPosition = Settings::GetDouble("ARCGeometry/Radius")
                         + Settings::GetDouble("RadiatorCell/CoolingThickness");
  if(CellRowNumber == 0 && CellColumnNumber == 0) {
    return Vector(XPosition, 0.0, 0.0);
  }
  if(CellRowNumber == 1) {
    const double ZPosition = HexagonSize*static_cast<double>(CellColumnNumber);
    return Vector(XPosition, 0.0, ZPosition);
  } else if(CellRowNumber == 2) {
    const double ZPosition = HexagonSize*(static_cast<double>(CellColumnNumber) - 0.5);
    const double BarrelRadius = Settings::GetDouble("ARCGeometry/Radius")
                              + Settings::GetDouble("RadiatorCell/CoolingThickness");
    const double DeltaPhi = 0.5*TMath::Sqrt(3)*HexagonSize/BarrelRadius;
    const ROOT::Math::RotationZ PhiRotation(DeltaPhi);
    return PhiRotation(Vector(XPosition, 0.0, ZPosition));
  } else {
    throw std::invalid_argument("Invalid cell row number: "
				+ std::to_string(CellRowNumber));
  }
}

Rotation3D BarrelRadiatorCell::GetCellOrientation(double HexagonSize,
						  std::size_t CellRowNumber) {
  const RotationY Rotation1(-TMath::Pi()/2.0);
  const RotationZ Rotation2(-TMath::Pi());
  if(CellRowNumber == 1) {
    return Rotation2*Rotation1;
  } else {
    const double BarrelRadius = Settings::GetDouble("ARCGeometry/Radius")
                              + Settings::GetDouble("RadiatorCell/CoolingThickness");
    const double DeltaPhi = 0.5*TMath::Sqrt(3)*HexagonSize/BarrelRadius;
    const ROOT::Math::RotationZ Rotation0(-DeltaPhi);
    return Rotation2*Rotation1*Rotation0;
  }
}
