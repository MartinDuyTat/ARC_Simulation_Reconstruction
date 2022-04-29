// Martin Duy Tat 29th April 2022
/**
 * The TrackingVolume is the inner detector that contains the magnetic field and the interaction point
 */

#ifndef TRACKINGVOLUME
#define TRACKINGVOLUME

class TrackingVolume {
 public:
  /**
   * Constructor that sets up the geometry of the tracking volume
   * @param Radius The inner radius of the tracker
   * @param FieldStrength The magnetic field strenth, in T
   */
  TrackingVolume(double Radius, double FieldStrength);
  /**
   * Get inner tracking radius
   */
  double GetRadius() const;
  /**
   * Get magnetic field strength
   */
  double GetFieldStrength() const;
 private:
  /**
   * The inner tracking radius
   */
  const double m_Radius;
  /**
   * The magnetic field strength
   */
  const double m_FieldStrength;
};

#endif
