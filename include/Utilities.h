// Martin Duy Tat 6th September 2022
/**
 * Various utility functions
 */

#ifndef UTILITIES
#define UTILITIES

#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"RadiatorCell.h"

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
   * Swap the x and z coordinates if the radiator cell is an end cap radiator cell
   * @param radiatorCell The radiator cell
   * @param v The vector we want to swap
   */
  Vector SwapXZForEndCap(const RadiatorCell *radiatorCell, Vector v);
};

#endif
