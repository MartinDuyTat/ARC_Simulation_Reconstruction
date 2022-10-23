// Martin Duy Tat 27th August 2022

#include"TArc.h"
#include"Math/RotationZ.h"
#include"TLine.h"
#include"HalfRadiatorCell.h"
#include"Settings.h"

using RotationZ = ROOT::Math::RotationZ;

HalfRadiatorCell::HalfRadiatorCell(std::size_t CellColumnNumber,
				   std::size_t CellRowNumber,
				   double HexagonSize):
  BarrelRadiatorCell(CellColumnNumber, CellRowNumber, HexagonSize) {
}

bool HalfRadiatorCell::IsInsideCell(const ARCVector &Position) const {
  if(Position.LocalVector().X() > 0.0) {
    return false;
  } else {
    return BarrelRadiatorCell::IsInsideCell(Position);
  }
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
HalfRadiatorCell::DrawRadiatorGeometry() const {
  // First find the mirror intersections with the walls in local coordinates
  // by solving a quadratic
  const auto MirrorCentre = GetMirrorCentre();
  const double b = MirrorCentre.Z();
  const double c_left = TMath::Power(m_HexagonSize/2.0, 2)
                      + MirrorCentre.Mag2()
                      - m_MirrorCurvature*m_MirrorCurvature
                      + MirrorCentre.X()*m_HexagonSize;
  const double c_right = MirrorCentre.Mag2()
                       - m_MirrorCurvature*m_MirrorCurvature;
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
  const double s2 = bb + TMath::Sqrt(bb*bb - cc);
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
  TLine RightLine(m_Position.Z(),
		  Settings::GetDouble("ARCGeometry/Radius") - m_VesselThickness,
		  m_Position.Z(),
		  Settings::GetDouble("ARCGeometry/Radius")
		  + m_RadiatorThickness - m_VesselThickness);
  RightLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(RightLine), ""));
  // Draw the detector plane
  TLine DetectorLine(-m_HexagonSize/2.0 + m_Position.Z(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_CoolingThickness,
		     m_Position.Z(),
		     Settings::GetDouble("ARCGeometry/Radius") + m_CoolingThickness);
  DetectorLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(DetectorLine), ""));
  // Draw the aerogel plane
  TLine AerogelLine(-m_HexagonSize/2.0 + m_Position.Z(),
		    Settings::GetDouble("ARCGeometry/Radius")
		    + m_CoolingThickness + m_AerogelThickness,
		    m_Position.Z(),
		    Settings::GetDouble("ARCGeometry/Radius")
		    + m_CoolingThickness + m_AerogelThickness);
  AerogelLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(AerogelLine), ""));
  // Draw the bottom cooling plane
  TLine CoolingLine(-m_HexagonSize/2.0 + m_Position.Z(),
		    Settings::GetDouble("ARCGeometry/Radius"),
		    m_Position.Z(),
		    Settings::GetDouble("ARCGeometry/Radius"));
  CoolingLine.SetLineColor(kBlack);
  Objects.push_back(std::make_pair(std::make_unique<TLine>(CoolingLine), ""));
  // Draw SiPM
  const auto RadiatorPosition = ReversePhi(GetRadiatorPosition());
  Objects.push_back(std::make_pair(m_Detector.DrawSiPM(RadiatorPosition), ""));
  return Objects;
}
