// Martin Duy Tat 23rd June 2022
/**
 * BarrelRadiatorArray contains the ARC barrel geometry of radiators
 */

#ifndef BARRELRADIATORARRAY
#define BARRELRADIATORARRAY

#include<vector>
#include<memory>
#include<utility>
#include<string>
#include"TObject.h"
#include"Math/Vector3Dfwd.h"
#include"RadiatorCell.h"
#include"RadiatorArray.h"

using Vector = ROOT::Math::XYZVector;

class ParticleTrack;

class BarrelRadiatorArray: public RadiatorArray {
 public:
  /**
   * Constructor that sets up the geometry
   */
  BarrelRadiatorArray();
  /**
   * Delete copy constructor
   */
  BarrelRadiatorArray(const BarrelRadiatorArray &radiatorCell) = delete;
  /**
   * Draw all the radiator cells
   */
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
  DrawRadiatorArray() const override;
  /**
   * Check which radiator the particle goes through
   * It will map the track momentum and position if the track hits an equivalent radiator cell
   */
  virtual const RadiatorCell* FindRadiator(ParticleTrack &particleTrack) const override;
  /**
   * Check if particle track hits below the main row
   * @param x Position along the z direction of the barrel, or x in local coordinates
   * @param y Position along the curved direction, or y in location coordinates
   */
  bool IsBelowMainRow(double x, double y) const;
  /**
   * Check if particle track hits above the upper row
   * @param x Position along the z direction of the barrel, or x in local coordinates
   * @param y Position along the curved direction, or y in location coordinates
   * @param Radius The radial coordinate of the particle in global coordinates
   */
  bool IsAboveUpperRow(double x, double y, double Radius) const;
  /**
   * Check if particle track hits above the main row
   * @param x Position along the z direction of the barrel, or x in local coordinates
   * @param y Position along the curved direction, or y in location coordinates
   */
  bool IsAboveMainRow(double x, double y) const;
 private:
  /**
   * Flag that is true if the full geometry is included, otherwise a single cell in the middle is used
   */
  const bool m_FullArray;
  /**
   * Radius of the barrel
   */
  const double m_BarrelRadius;
  /**
   * The azimuthal rotation between each two rows
   */
  const double m_DeltaPhi;
  /**
   * Helper function to find the index for the radiator cell
   */
  virtual int FindRadiatorIndex(std::size_t i, std::size_t j) const override;
};

#endif
