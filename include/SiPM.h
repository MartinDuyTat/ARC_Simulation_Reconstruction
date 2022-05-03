// Martin Duy Tat 29th April 2022
/**
 * SiPM represents a single Silicon Photo-Multiplier, which is located on the detector plane of each radiator cell
 * All positions are in the local radiator coordinate system
 */

#ifndef SIPM
#define SIPM

#include<vector>
#include<string>
#include"Photon.h"

struct PhotonHit {
  /**
   * Constructor saving the photon hits
   */
  PhotonHit(double x, double y, const Photon *photon): x(x), y(y), m_Photon(photon) {}
  /**
   * Detector hit coordinates
   */
  const double x;
  const double y;
  /**
   * Pointer to the photon that caused this hit
   */
  const Photon *m_Photon;
};

class SiPM {
 public:
  /**
   * Constructor setting up the detector coordinates and detector size
   */
  SiPM();
  /**
   * Add a photon hit
   */
  void AddPhotonHit(const Photon &photon);
  /**
   * Plot photon hits
   */
  void PlotHits(const std::string &Filename) const;
  /**
   * Get photon hits
   */
  const std::vector<PhotonHit>& GetPhotonHits() const;
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

#endif