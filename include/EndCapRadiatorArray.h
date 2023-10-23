// Martin Duy Tat 23rd June 2022
/**
 * EndCapRadiatorArray contains the ARC end cap geometry of radiators
 */

#ifndef ENDCAPRADIATORARRAY
#define ENDCAPRADIATORARRAY

#include<vector>
#include<memory>
#include<utility>
#include<string>
#include"TObject.h"
#include"Math/Vector3Dfwd.h"
#include"RadiatorCell.h"
#include"RadiatorArray.h"

using Vector = ROOT::Math::XYZVector;

class Particle;

class EndCapRadiatorArray: public RadiatorArray {
 public:
  /**
   * Constructor that sets up the geometry
   */
  EndCapRadiatorArray();
  /**
   * Delete copy constructor
   */
  EndCapRadiatorArray(const EndCapRadiatorArray &radiatorCell) = delete;
  /**
   * Draw all the radiator cells
   */
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
  DrawRadiatorArray() const override;
  /**
   * Check which radiator the particle goes through
   * It will map the track momentum and position if the track hits an equivalent radiator cell
   */
  virtual const RadiatorCell* FindRadiator(Particle &particle) const override;
  /**
   * Check if particle track hits below a particular row
   * @param Row The row number
   * @param x Position along the z direction of the barrel, or x in local coordinates
   * @param y Position along the curved direction, or y in location coordinates
   */
  bool IsBelowThisRow(std::size_t Row, double x, double y) const;
  /**
   * Check if particle track hits below a row with odd row number
   * It is assumed the y cooridnate has been shifted so it's inside the main row
   * @param x Position along the z direction of the barrel, or x in local coordinates
   * @param y Position along the curved direction, or y in location coordinates
   */
  bool IsBelowOddRow(double x, double y) const;
  /**
   * Check if particle track hits below a row with even row number
   * It is assumed the y cooridnate has been shifted so it's inside the main row
   * @param x Position along the z direction of the barrel, or x in local coordinates
   * @param y Position along the curved direction, or y in location coordinates
   */
  bool IsBelowEvenRow(double x, double y) const;
 private:
  /**
   * Inner radius of end cap
   */
  const double m_InnerRadius;
  /**
   * Outer radius of end cap
   */
  const double m_OuterRadius;
  /**
   * Helper function to find the index for the radiator cell
   */
  virtual int FindRadiatorIndex(std::size_t i, std::size_t j) const override;
};

#endif
