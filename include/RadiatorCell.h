// Martin Duy Tat 29th April 2022
/**
 * RadiatorCell describes the geometry of a single radiator cell of the ARC
 * The radiator cell has a local coordinate system with its origin in the middle of the detector plane
 * All lengths are in meter
 */

#ifndef RADIATORCELL
#define RADIATORCELL

#include<string>
#include<vector>
#include<utility>
#include<memory>
#include"TMath.h"
#include"Math/Vector3Dfwd.h"
#include"SiPM.h"
#include"Photon.h"

using Vector = ROOT::Math::XYZVector;

class RadiatorCell {
 public:
  /**
   * Constructor that sets up the geometry with a default cell position
   * @param CellColumnNumber The column number of this cell
   * @param CellRowNumber The row number of this cell
   * @param HexagonSize The length between two opposide edges (not points) of the hexagons
   * @param Prefix Prefix in the names of cells in options file
   */
  RadiatorCell(int CellColumnNumber,
	       int CellRowNumber,
	       double HexagonSize,
	       const std::string &Prefix = "");
  /**
   * Need virtual destructor because of polymorphism
   */
  virtual ~RadiatorCell() = default;
  /**
   * Get total radiator cell thickness
   */
  double GetRadiatorThickness() const;
  /**
   * Get thickness of vessel
   */
  double GetVesselThickness() const;
  /**
   * Get thickness of cooling plate
   */
  double GetCoolingThickness() const;
  /**
   * Get thickness of aerogel
   */
  double GetAerogelThickness() const;
  /**
   * Get the mirror centre of curvature
   */
  const Vector& GetMirrorCentre() const;
  /**
   * Get the mirror curvature
   */
  double GetMirrorCurvature() const;
  /**
   * Get the radiator cell position in the global coordinates
   */
  virtual const Vector& GetRadiatorPosition() const = 0;
  /**
   * Function that checks if the position is inside the radiator
   */
  virtual bool IsInsideCell(const Vector &Position) const = 0;
  /**
   * Function that checks if the photon is inside the radiator
   */
  bool IsInsideCell(const Photon &photon) const;
  /**
   * Draw radiator geometry
   */
  virtual std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
  DrawRadiatorGeometry() const = 0;
  /**
   * Get length between two edges of a hexagon cell
   */
  double GetHexagonSize() const;
  /**
   * Get cell number
   */
  std::pair<int, int> GetCellNumber() const;
  /**
   * Get the detector
   */
  const SiPM& GetDetector() const;
  /**
   * Set the mirror curvature
   * @param Curvature New mirror curvature
   */
  void SetMirrorCurvature(double Curvature);
  /**
   * Set the mirror x position, relative to the default position
   */
  virtual void SetMirrorXPosition(double x);
  /**
   * Set the mirror z position, relative to the default position
   */
  void SetMirrorZPosition(double z);
  /**
   * Set the detector X position
   */
  void SetDetectorPosition(double x);
  /**
   * Set the detector tilt, in radians
   */
  void SetDetectorTilt(double Angle);
 protected:
  /**
   * Get the centre of curvature of the mirror z coordinate in local coordinates
   */
  double GetMirrorCurvatureCentreZ() const;
  /**
   * The full thickness of the radiator cell
   */
  const double m_RadiatorThickness;
  /**
   * Thickness of the vessel
   */
  const double m_VesselThickness;
  /**
   * Thickness of cooling plate
   */
  const double m_CoolingThickness;
  /**
   * Thickness of the aerogel
   */
  const double m_AerogelThickness;
  /**
   * The length between two edges of a hexagon
   */
  const double m_HexagonSize;
  /**
   * Radiator cell row number (first) and column number (second)
   * For a single central cell, it is assigned number (0, 0)
   * For an array of radiator cells, numbering starts from (0, 1) from the middle of the main row
   */
  const std::pair<int, int> m_CellNumber;
  /**
   * The SiPM in the radiator cell
   */
  SiPM m_Detector;
  /**
   * Radius of curvature of the spherical mirror
   */
  double m_MirrorCurvature;
  /**
   * Default centre of curvature of the mirror
   */
  const Vector m_DefaultMirrorCentre;
  /**
   * Centre of curvature of the mirror
   */
  Vector m_MirrorCentre;
};

#endif
