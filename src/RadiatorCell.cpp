// Martin Duy Tat 29th April 2022

#include<vector>
#include<utility>
#include<string>
#include<memory>
#include<stdexcept>
#include"Math/Vector3Dfwd.h"
#include"Math/Rotation3D.h"
#include"RadiatorCell.h"
#include"Settings.h"
#include"Photon.h"

RadiatorCell::RadiatorCell(std::size_t CellColumnNumber,
			   std::size_t CellRowNumber,
			   double HexagonSize,
			   const Vector &Position,
			   const Rotation3D &Rotation,
			   const std::string &Prefix):
  m_Position(Position),
  m_Rotation(Rotation),
  m_RadiatorThickness(Settings::GetDouble("RadiatorCell/RadiatorThickness")),
  m_VesselThickness(Settings::GetDouble("RadiatorCell/VesselThickness")),
  m_CoolingThickness(Settings::GetDouble("RadiatorCell/CoolingThickness")),
  m_AerogelThickness(Settings::GetDouble("RadiatorCell/AerogelThickness")),
  m_HexagonSize(HexagonSize),
  m_CellNumber(std::make_pair(CellColumnNumber, CellRowNumber)),
  m_MirrorCurvature(Settings::GetDouble("RadiatorCell/MirrorCurvature")),
  m_DefaultMirrorCentre(0.0, 0.0, GetMirrorCurvatureCentreZ()),
  m_MirrorCentre(m_DefaultMirrorCentre, m_Position, m_Rotation) {
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
  return m_MirrorCentre.LocalVector();
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
  return IsInsideCell(photon.GetPosition());
}

const Vector& RadiatorCell::GetRadiatorPosition() const {
  return m_Position;
}

const Rotation3D& RadiatorCell::GetRadiatorRotation() const {
  return m_Rotation;
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

bool RadiatorCell::IsEdgeCell() const {
  return (m_CellNumber == std::make_pair(std::size_t{8}, std::size_t{1}) ||
         m_CellNumber == std::make_pair(std::size_t{9}, std::size_t{2}));
}

RotationZ RadiatorCell::ReversePhiRotation() const {
  return RotationZ(0.0);
}
