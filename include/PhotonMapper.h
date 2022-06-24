// Martin Duy Tat 1st May 2022
/**
 * PhotonMapper is a namespace containing the functions for mapping photons from their emission point to the detector plane
 */

#ifndef PHOTONMAPPER
#define PHOTONMAPPER

#include"RadiatorCell.h"
#include"Photon.h"

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
  void TracePhotonToDetector(Photon &photon);
  /**
   * Trace photon from emission point to detector and register detector hit
   */
  void TracePhoton(Photon &photon);
}

#endif
