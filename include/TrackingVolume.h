// Martin Duy Tat 29th April 2022
/**
 * The TrackingVolume is the inner detector that contains the magnetic field and the interaction point
 */

#ifndef TRACKINGVOLUME
#define TRACKINGVOLUME

#include<memory>
#include<utility>
#include<string>
#include<vector>
#include"TObject.h"

class TrackingVolume {
 public:
  /**
   * Constructor that sets up the geometry of the tracking volume
   */
  TrackingVolume();
  /**
   * Get inner tracking radius
   */
  double GetRadius() const;
  /**
   * Get magnetic field strength
   */
  double GetFieldStrength() const;
  /**
   * Draw ARC geometry
   */
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> DrawARCGeometry() const;
 private:
  /**
   * The inner tracking radius
   */
  const double m_Radius;
  /**
   * Length of detector
   */
  const double m_Length;
  /**
   * Number of cells in theta direction
   */
  const int m_ThetaCells;
  /**
   * Number of cells in phi direction
   */
  const int m_PhiCells;
  /**
   * The magnetic field strength
   */
  const double m_FieldStrength;
};

#endif
