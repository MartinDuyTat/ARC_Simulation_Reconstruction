// Martin Duy Tat 29th April 2022
/**
 * SiPM represents a single Silicon Photo-Multiplier, which is located on the detector plane of each radiator cell
 * All positions are in the local radiator coordinate system
 */

#ifndef SIPM
#define SIPM

struct PhotonHit {
  /**
   * Constructor saving the photon hits
   */
  PhotonHit(double x, double y): x(x), y(y) {}
  /**
   * Detector hit coordinates
   */
  const double x;
  const double y;
};

class SiPM {
 public:
  /**
   * Constructor setting up the detector coordinates and detector size
   */
  SiPM();
 private:
  /**
   * Size of detector in x direction
   */
  const double m_DetectorSizeX;
  /**
   * Size of detector in y direction
   */
  const double m_DetectorSizeY;
  /**
   * Centre of detector in x direction
   */
  const double m_DetectorPositionX;
  /**
   * Centre of detector in y direction
   */
  const double m_DetectorPositionY;
  /**
   * Vector of photon hits
   */
  std::vector<PhotonHit> m_PhotonHits;
};
