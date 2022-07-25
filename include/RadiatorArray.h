// Martin Duy Tat 23rd June 2022
/**
 * RadiatorArray contains the full ARC geometry of radiators
 */

#ifndef RADIATORARRAY
#define RADIATORARRAY

#include<vector>
#include<memory>
#include<utility>
#include<string>
#include"TObject.h"
#include"Math/Vector3Dfwd.h"
#include"RadiatorCell.h"

using Vector = ROOT::Math::XYZVector;

class ParticleTrack;

class RadiatorArray {
 public:
  /**
   * Constructor that sets up the geometry
   */
  RadiatorArray();
  /**
   * Draw all the radiator cells
   */
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> DrawRadiatorArray() const;
  /**
   * Operator overload to get the individual radiator cells
   */
  RadiatorIter operator()(int i, int j);
  /**
   * Check which radiator the particle goes through
   * It will map the track momentum and position if the track hits an equivalent radiator cell
   */
  RadiatorIter FindRadiator(ParticleTrack &particleTrack);
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
   */
  bool IsAboveUpperRow(double x, double y) const;
  /**
   * Check if particle track hits above the main row
   * @param x Position along the z direction of the barrel, or x in local coordinates
   * @param y Position along the curved direction, or y in location coordinates
   */
  bool IsAboveMainRow(double x, double y) const;
  /**
   * Reset the detectors of all radiator cells
   */
  void ResetDetectors();
 private:
  /**
   * Vector containing all the radiator cells
   */
  std::vector<RadiatorCell> m_Cells;
  /**
   * Flag that is true if the full geometry is included, otherwise a single cell in the middle is used
   */
  const bool m_FullArray;
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
