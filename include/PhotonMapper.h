// Martin Duy Tat 1st May 2022
/**
 * PhotonMapper is a namespace containing the functions for mapping photons from their emission point to the detector plane
 */

#ifndef PHOTONMAPPER
#define PHOTONMAPPER

#include<memory>
#include"RadiatorCell.h"
#include"Photon.h"
#include"SiPM.h"

namespace PhotonMapper {
  /**
   * Find the distance between the photon and the mirror
   */
  double PhotonMirrorDistance(const Photon &photon);
  /**
   * Trace photon from emission point to mirror
   */
  void TracePhotonToMirror(Photon &photon);
  /**
   * Trace photon from mirror to detector plane and register detector hit
   */
  PhotonHit TracePhotonToDetector(Photon &photon);
  /**
   * Trace photon from emission point to detector and register detector hit
   */
  std::unique_ptr<PhotonHit> TracePhoton(Photon &photon);
  /**
   * Check if photon has been scattered in the aerogel
   */
  bool IsScatteredInAerogel(const Photon &photon);
}

#endif
