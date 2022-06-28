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
using RadiatorIter = std::vector<RadiatorCell>::iterator;

class RadiatorCell {
 public:
  /**
   * enum class describing if the radiator cell is the first, middle or last cell
   */
  enum class FirstMiddleLast{First, Middle, Last, Single};
  /**
   * Constructor that sets up the geometry
   * @param CellNumber Unique number that identifies the position of the cell
   */
  RadiatorCell(int CellNumber);
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
  const Vector& GetRadiatorPosition() const;
  /**
   * Function that checks if the position is inside the radiator in the theta direction
   */
  bool IsInsideThetaBoundary(const Vector &Position) const;
  /**
   * Function that checks if the photon is inside the radiator in the theta direction
   */
  bool IsInsideThetaBoundary(const Photon &photon) const;
  /**
   * Function that checks if the photon is inside the radiator in the phi direction
   */
  bool IsInsidePhiBoundary(const Photon &photon) const;
  /**
   * Draw radiator geometry
   */
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> DrawRadiatorGeometry() const;
  /**
   * Get length of cell in theta direction
   */
  double GetThetaLength() const;
  /**
   * Get cell number
   */
  double GetCellNumber() const;
  /**
   * Get the radiator cell position to be first, middle or last
   */
  FirstMiddleLast GetFirstMiddleLast() const;
  /**
   * Get the detector
   */
  SiPM& GetDetector();
 private:
  /**
   * Get the centre of curvature of the mirror z coordinate in local coordinates
   */
  double GetMirrorCurvatureCentreZ() const;
  /**
   * Length of the radiator cell in the theta direction
   */
  const double m_ThetaLength;
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
   * Position of radiator cell in global coordinates, but rotated so that x is along the theta direction and y is along the phi direction
   */
  const Vector m_Position;
  /**
   * The SiPM in the radiator cell
   */
  SiPM m_Detector;
  /**
   * Radius of curvature of the spherical mirror
   */
  const double m_MirrorCurvature;
  /**
   * Centre of curvature of the mirror
   */
  const Vector m_MirrorCentre;
  /**
   * Angular size of the radiator cell in the phi direction
   */
  const double m_DeltaPhi;
  /**
   * Radiator cell number
   * For a single central cell, it is assigned number 0
   * For an array of radiator cells, numbering starts from 1 from the middle to the right and -1 to the left
   */
  const int m_CellNumber;
  /**
   * Helper function to get cell position based on cell number
   */
  Vector GetCellPosition(int CellNumber) const;
  /**
   * Label that describes if the radiator cell is the first, last or if it's one of the middle radiators
   */
  FirstMiddleLast m_FirstMiddleLast;
  /**
   * Helper function to determine the mirror curvature
   */
  double DetermineMirrorCurvature() const;
  /**
   * Helper function to figure out position of SiPM
   */
  double DetermineSiPMPositionX() const;
};

#endif
