// Martin Duy Tat 3rd May 2022
/**
 * ReconstructedPhoton contains the reconstructed Cherenkov angle from a photon hit
 */

#ifndef RECONSTRUCTEDPHOTON
#define RECONSTRUCTEDPHOTON

#include"Photon.h"

struct ReconstructedPhoton {
  /**
   * Trivial constructor
   * @param photon The original photon
   */
  ReconstructedPhoton(const Photon &photon);
  /**
   * The photon that is being reconstructed
   */
  const Photon *m_Photon;
  /**
   * The reconstructed Cherenkov angle, using the true emission point
   */
  double m_CosCherenkovAngle_TrueEmissionPoint;
  /**
   * The reconstructed Cherenkov angle, assming emission in the middle
   */
  double m_CosCherenkovAngle;
};

#endif
