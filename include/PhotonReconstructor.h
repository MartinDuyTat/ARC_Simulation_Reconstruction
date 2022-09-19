// Martin Duy Tat 3rd May 2022
/**
 * PhotonReconstructer is a namespace with all the tools for reconstructing photons by solving a quartic equation
 */

#ifndef PHOTONRECONSTRUCTOR
#define PHOTONRECONSTRUCTOR

#include<array>
#include"Math/DisplacementVector3D.h"
#include"Photon.h"
#include"ParticleTrack.h"
#include"RadiatorCell.h"
#include"ReconstructedPhoton.h"

using Vector = ROOT::Math::XYZVector;

namespace PhotonReconstructor {

  /**
   * Reconstruct cosine of Cherenkov angle of photon from the photon hit and charged track
   * @param Particle The charged particle that emitted the photon
   * @param photonHit The photon hit in the detector
   * @param TrueEmissionPoint Set to true to use the true emission point
   * @param Radiator The radiator that the photon was emitted from
   */
  double ReconstructCosCherenkovAngle(const ParticleTrack &Particle,
				      const PhotonHit &photonHit,
				      bool TrueEmissionPoint,
				      Photon::Radiator Radiator);
  /**
   * Reconstruct photon from the photon hit, charged track and mirror geometry
   * @param Particle The charged particle that emitted the photon
   * @param photonHit The photon hit in the detector
   * @param radiatorCell The radiator cell
   * @param Radiator The radiator that the photon was emitted from
   */
  ReconstructedPhoton ReconstructPhoton(const ParticleTrack &Particle,
					const PhotonHit &photonHit,
					Photon::Radiator Radiator);
  /**
   * Struct with solutions to the quartic equation
   */
  struct QuarticSolution {
    /**
     * Sine and cosine of angle between emission and reflection point, relative to the mirror centre
     */
    std::array<double, 2> m_SinBeta;
    std::array<double, 2> m_CosBeta;
    /**
     * Flag that is true if there are more than 2 solutions
     */
    bool m_DegenerateSolution = false;
  };
  /**
   * Solve quartic equation to determine sine of angle between emission and reflection point, relative to the mirror centre
   * @param EmissionMirrorDist Distance between the emission point and mirror centre
   * @param DetectionMirrorParaDist Parallel component of distance between detection point and mirror centre
   * @param DetectionMirrorPerpDist Perpendicular component of distance between detection point and mirror centre
   */
  QuarticSolution SolveQuartic(double EmissionMirrorDist,
			       double DetectionMirrorParaDist,
			       double DetectionMirrorPerpDist);
  /**
   * Helper function to calculate the reflection point after solving the quartic equation
   */
  Vector GetReflectionPoint(const Vector &EmissionPoint,
			    const Vector &DetectionPoint,
			    double DetectionMirrorParaDist,
			    double DetectionMirrorPerpDist,
			    double EmissionMirrorDist,
			    const QuarticSolution &quarticSolution,
			    std::size_t SolutionNumber);
}

#endif
