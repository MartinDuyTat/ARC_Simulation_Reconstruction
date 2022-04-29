// Martin Duy Tat 29th April 2022
/**
 * RadiatorCell describes the geometry of a single radiator cell of the ARC
 * The radiator cell has a local coordinate system with its origin in the middle of the detector plane
 * All lengths are in meter
 */

#ifndef RADIATORCELL
#define RADIATORCELL

#include"Vector3Dfwd.h"

using Vector = ROOT::Math::XYZVector;

class RadiatorCell {
 public:
  /**
   * Constructor that sets up the geometry
   * @param Position Position of radiator cell in global coordinates
   */
  RadiatorCell(const Vector &Position);
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
 private:
  /**
   * Get the centre of curvature of the mirror z coordinate in local coordinates
   */
  double GetMirrorCurvatureCentreZ() const;
  /**
   * Get the radiator cell position in the global coordinates
   */
  const Vector& GetRadiatorPosition() const;
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
   * Radius of curvature of the spherical mirror
   */
  const double m_MirrorCurvature;
  /**
   * Centre of curvature of the mirror
   */
  const Vector m_MirrorCentre;
  /**
   * Position of origin of radiator cell in global coordinate system
   */
  const Vector m_Position;
};

#endif
