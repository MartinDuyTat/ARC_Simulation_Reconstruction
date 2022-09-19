// Martin Duy Tat 29th April 2022

#include<vector>
#include<utility>
#include<string>
#include<memory>
#include<stdexcept>
#include"RadiatorCell.h"
#include"Settings.h"
#include"Photon.h"

RadiatorCell::RadiatorCell(std::size_t CellColumnNumber,
			   std::size_t CellRowNumber,
			   double HexagonSize,
			   const std::string &Prefix):
  m_RadiatorThickness(Settings::GetDouble("RadiatorCell/RadiatorThickness")),
  m_VesselThickness(Settings::GetDouble("RadiatorCell/VesselThickness")),
  m_CoolingThickness(Settings::GetDouble("RadiatorCell/CoolingThickness")),
  m_AerogelThickness(Settings::GetDouble("RadiatorCell/AerogelThickness")),
  m_HexagonSize(HexagonSize),
  m_CellNumber(std::make_pair(CellColumnNumber, CellRowNumber)),
  m_Detector(),
  m_MirrorCurvature(Settings::GetDouble("RadiatorCell/MirrorCurvature")),
  m_DefaultMirrorCentre(0.0, 0.0, GetMirrorCurvatureCentreZ()),
  m_MirrorCentre(m_DefaultMirrorCentre) {
  const std::string RadiatorName = Prefix + "Radiator_c"
                                 + std::to_string(m_CellNumber.first)
                                 + "_r"
                                 + std::to_string(m_CellNumber.second)
                                 + "_";
  const std::string CurvatureName = "RadiatorCell/" + RadiatorName + "Curvature";
  if(Settings::Exists(CurvatureName)) {
    SetMirrorCurvature(Settings::GetDouble(CurvatureName));
  }
  const std::string XPositionName = "RadiatorCell/" + RadiatorName + "XPosition";
  if(Settings::Exists(XPositionName)) {
    SetMirrorXPosition(Settings::GetDouble(XPositionName));
  }
  const std::string ZPositionName = "RadiatorCell/" + RadiatorName + "ZPosition";
  if(Settings::Exists(ZPositionName)) {
    SetMirrorZPosition(Settings::GetDouble(ZPositionName));
  }
  const std::string DetPositionName = "RadiatorCell/" + RadiatorName + "DetPosition";
  if(Settings::Exists(DetPositionName)) {
    SetDetectorPosition(Settings::GetDouble(DetPositionName));
  }
  const std::string DetTiltName = "RadiatorCell/" + RadiatorName + "DetTilt";
  if(Settings::Exists(DetTiltName)) {
    SetDetectorTilt(Settings::GetDouble(DetTiltName));
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
  return m_RadiatorThickness
       - m_MirrorCurvature
       - 2*m_VesselThickness
       - m_CoolingThickness;
}

double RadiatorCell::GetMirrorCurvature() const {
  return m_MirrorCurvature;
}

bool RadiatorCell::IsInsideCell(const Photon &photon) const {
  return IsInsideCell(photon.m_Position);
}

double RadiatorCell::GetHexagonSize() const {
  return m_HexagonSize;
}

std::pair<std::size_t, std::size_t> RadiatorCell::GetCellNumber() const {
  return m_CellNumber;
}

const SiPM& RadiatorCell::GetDetector() const {
  return m_Detector;
}

void RadiatorCell::SetMirrorCurvature(double Curvature) {
  m_MirrorCurvature = Curvature;
}

void RadiatorCell::SetMirrorXPosition(double x) {
  m_MirrorCentre.SetX(m_DefaultMirrorCentre.X() + x);
}

void RadiatorCell::SetMirrorZPosition(double z) {
  m_MirrorCentre.SetZ(m_DefaultMirrorCentre.Z() + z);
}

void RadiatorCell::SetDetectorPosition(double x) {
  m_Detector.SetDetectorPosition(x);
}

void RadiatorCell::SetDetectorTilt(double Angle) {
  m_Detector.SetDetectorTilt(Angle);
}

bool RadiatorCell::IsDetectorInsideCell() const {
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
