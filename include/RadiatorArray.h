// Martin Duy Tat 2nd September 2022
/**
 * RadiatorArray is an abstract class that represent an array of radiator cells
 */

#ifndef RADIATORARRAY
#define RADIATORARRAY

#include<vector>
#include<memory>
#include"RadiatorCell.h"

class ParticleTrack;

class RadiatorArray {
 public:
  /**
   * Default constructor
   */
  RadiatorArray();
  /**
   * Need virtual destructor because of polymorphism
   */
  virtual ~RadiatorArray() = default;
  /**
   * Operator overload to get the individual radiator cells
   */
  virtual const RadiatorCell* operator()(int i, int j) = 0;
  /**
   * Check which radiator the particle goes through
   * It will map the track momentum and position if the track hits an equivalent radiator cell
   */
  virtual const RadiatorCell* FindRadiator(ParticleTrack &particleTrack) = 0;
  /**
   * Draw all the radiator cells
   */
  virtual std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
  DrawRadiatorArray() const = 0;
 protected:
  /**
   * Vector containing all the radiator cells
   */
  std::vector<std::unique_ptr<RadiatorCell>> m_Cells;
  /**
   * Number of cells in theta direction along the main row
   */
  const int m_NumberMainRowCells;
  /**
   * The horizontal distance between two hexagons
   */
  const double m_xHexDist;
  /**
   * The vertical distance between two hexagons
   */
  const double m_yHexDist;
};

#endif
