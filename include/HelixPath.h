// Martin Duy Tat 28th October 2022
/**
 * HelixPath is a class that stores the helical trajectory of a charged particle
 * in a magnetic field
 */

#ifndef HELIXPATH
#define HELIXPATH

#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"HelixFunctor.h"

using Vector = ROOT::Math::XYZVector;
using RotationZ = ROOT::Math::RotationZ;

class HelixPath {
 public:
  /**
   * Constructor that takes in the particle charge, momentum and magnetic field
   * @param Momentum Initial momentum vector of charged particle
   * @param Q Particle charge, +1 or -1
   * @param B Magnetic field
   */
  HelixPath(const Vector &Momentum, int Q, double B);
  /**
   * We don't need copy constructor
   */
  HelixPath(const HelixPath &Helix) = default;
  /**
   * Get the particle position from the path length
   * @param s Path length coordinate
   */
  Vector GetPosition(double s) const;
  /**
   * Get the particle direction unit vector from the path length
   * @param s Path length coordinate
   */
  Vector GetDirection(double s) const;
  /**
   * Find the path length that solves the equation f(x) = 0
   * @param Functor Function that we want to find root of
   */
  double SolvePathLength(const HelixFunctor &Functor,
			 double Min, double Max) const;
  /**
   * Function that resets the origin
   * @param NewOrigin_s The new path length origin
   */
  void ResetOrigin(double NewOrigin_s);
  /**
   * Rotate around the phi direction to map particle to a valid radiator cell
   * Warning! Need to call ResetOrigin first!
   * @param DeltaPhi Angle that is rotated in phi
   */
  void MapPhi(double DeltaPhi);
  /**
   * Reflect everything in the z-direction
   */
  void ReflectZ();
  /**
   * Reflect everything in the y-direction
   * Warning! Need to call ResetOrigin first!
   */
  void ReflectY();
 private:
  /**
   * Unit vector in the initial particle direction
   */
  Vector m_n;
  /**
   * Signed radius of curvature
   */
  const double m_Radius;
  /**
   * Magnitude of inverse of radius of curvature
   */
  const double m_InvRadius;
  /**
   * Flag that is true when there is a non-zero magnetic field
   */
  const bool m_BFieldOn;
  /**
   * The path length at the current origin
   */
  double m_Origin_s;
  /**
   * The position of the origin
   */
  Vector m_Origin;
};

#endif
