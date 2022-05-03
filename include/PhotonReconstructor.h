// Martin Duy Tat 3rd May 2022
/**
 * PhotonReconstructer is a namespace with all the tools for reconstructing photons by solving a quartic equation
 */

#ifndef PHOTONRECONSTRUCTOR
#define PHOTONRECONSTRUCTOR

#include<array>
#include"Photon.h"
#include"ParticleTrack.h"
#include"RadiatorCell.h"
#include"ReconstructedPhoton.h"

namespace PhotonReconstructor {

  /**
   * Reconstruct photon from the photon hit, charged track and mirror geometry
   * @param Particle The charged particle that emitted the photon
   * @param photonHit The photon hit in the detector
   * @param radiatorCell The radiator cell
   */
  ReconstructedPhoton ReconstructPhoton(const ParticleTrack &Particle, const PhotonHit &photonHit, const RadiatorCell &radiatorCell);
  /**
   * Struct with solutions to the quartic equation
   */
  struct QuarticSolution {
    /**
     * Sine and cosine of angle between emission and reflection point, relative to the mirror centre
     */
    std::array<double, 2> m_SinBeta;
    std::array<double, 2> m_CosBeta;
  };
  /**
   * Solve quartic equation to determine sine of angle between emission and reflection point, relative to the mirror centre
   * @param EmissionMirrorDist Distance between the emission point and mirror centre
   * @param DetectionMirrorParaDist Parallel component of distance between detection point and mirror centre
   * @param DetectionMirrorPerpDist Perpendicular component of distance between detection point and mirror centre
   */
  QuarticSolution SolveQuartic(double EmissionMirrorDist, double DetectionMirrorParaDist, double DetectionMirrorPerpDist);

}

#endif
