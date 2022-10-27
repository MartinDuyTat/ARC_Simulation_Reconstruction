// Martin Duy Tat 6th September 2022
/**
 * Various utility functions
 */

#ifndef UTILITIES
#define UTILITIES

#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"RadiatorCell.h"
#include"RadiatorArray.h"
#include"ParticleTrack.h"
#include"Photon.h"

using Vector = ROOT::Math::XYZVector;

namespace Utilities {
  /**
   * Create a vector from spherical coordinates
   */
  Vector VectorFromSpherical(double R, double CosTheta, double Phi);
  /**
   * Generate a random barrel track which is uniform in z and phi
   */
  Vector GenerateRandomBarrelTrack(double &CosTheta, double &Phi);
  /**
   * Generate a random barrel track which is uniform in z and phi
   */
  Vector GenerateRandomBarrelTrackZRange(double z_min, double z_max);
  /**
   * Generate a random end cap track which is uniform in x and y
   */
  Vector GenerateRandomEndCapTrack();
  /**
   * Get the magnitude of the particle momentum
   */
  double GetMomentumMag();
  /**
   * Struct containing sum of Cherenkov angle, sum of square of Cherenkov
   * angle and total number of photons
   */
  struct ResolutionStruct {
    /**
     * The Cherenkov angle
     */
    double x = 0.0;
    /**
     * The number of photons
     */
    int N = 0;
    /**
     * Set to true if particle hits the top wall of the radiator cell
     */
    bool HitTopWall = false;
    /**
     * Set to true if the particle hits the correct cell
     */
    bool HitCorrectCell = false;
    /**
     * The average photon hit position
     */
    Vector CentreHitDistance{0.0, 0.0, 0.0};
  };
  /**
   * Track photons through radiator cells and return the information about resolutions
   */
  ResolutionStruct TrackPhotons(ParticleTrack particleTrack,
				const RadiatorCell &radiatorCell,
				const RadiatorArray &radiatorArray);
  /**
   * Get index of fraction
   * @param Radiator Gas or aerogel
   * @param Energy The photon energy
   * @param 
   */
  double GetIndexRefraction(Photon::Radiator Radiator,
			    bool ChromaticDispersion,
			    double Energy = 0.0);
  /**
   * Get the predicted Cherenkov angle
   * @param Momentum The particle momentum
   * @param ID The particle ID
   * @param Radiator Gas or aerogel
   */
  double GetPredictedCherenkovAngle(double Momentum,
				    int ID,
				    Photon::Radiator Radiator);
  /**
   * Get the difference in Cherenkov angle between two mass hypotheses
   * @param Momentum The particle momentum
   * @param ID1 The first particle hypothesis
   * @param ID2 The second particle hypothesis
   * @param Radiator Gas or aerogel
   */
  double GetCherenkovAngleDifference(double Momentum,
				     int ID1,
				     int ID2,
				     Photon::Radiator Radiator);
}

#endif
