// Martin Duy Tat 28th October 2022
/**
 * HelixFunctor is an abstract class for functors to solve the particle position
 */

#ifndef HELIXFUNCTOR
#define HELIXFUNCTOR

#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"

using Vector = ROOT::Math::XYZVector;

class HelixFunctor {
 public:
  /**
   * Constructor is default
   */
  HelixFunctor() = default;
  /**
   * Virtual () operator
   * @param Vect Vector position
   */
  virtual double operator()(const Vector &Vect) const = 0;
};

/**
 * BarrelHelixFunctor is a functor for tracking the particle to a certain barrel radius
 */
class BarrelHelixFunctor: public HelixFunctor {
 public:
  /**
   * Constructor that sets up the radius that we want to track particle to
   * @param Radius The radius we're tracking the particle to
   */
  BarrelHelixFunctor(double Radius);
  /**
   * We do not need copy constructor
   */
  BarrelHelixFunctor(const BarrelHelixFunctor &Functor) = delete;
  /**
   * Virtual () operator
   * @param Vect Vector position
   */
  virtual double operator()(const Vector &Vect) const override;
 private:
  /**
   * The radius we're tracking the particle to
   */
  const double m_Radius;
};

/**
 * ZPlaneHelixFunctor is a functor for tracking the particle through various
 * the different layers of a barrel cell
 */

class ZPlaneHelixFunctor: public HelixFunctor {
 public:
  /**
   * Constructor that takes in the position of the cell and the z-coordinate of
   * the local coordinate that we want to track particle to
   * @param z The z-coordinate
   * @param CellPosition The position of the radiator cell
   * @param Barrel Flag that is true for barrel cells, false for end cap cells
   */
  ZPlaneHelixFunctor(double z, const Vector &CellPosition, bool Barrel);
  /**
   * We do not need copy constructor
   */
  ZPlaneHelixFunctor(const ZPlaneHelixFunctor &Functor) = delete;
  /**
   * Virtual () operator
   * @param Vect Vector position
   */
  virtual double operator()(const Vector &Vect) const override;
 private:
  /**
   * The z-coordinate we're tracking the particle to
   */
  const double m_z;
  /**
   * The unit vector along the z-direction in local coordinates
   */
  const Vector m_UnitZ;
  /**
   * The z-distance from the detector centre to the radiator cell origin
   */
  const double m_CellPositionZ;
  /**
   * Helper function that determines the unit vector along the local z-coordinate
   * @param Barrel Flag that is true for barrel cells, false for end cap cells
   * @param Position Position of cell
   */
  static Vector GetUnitZ(bool Barrel, const Vector &Position);
};

/**
 * MirrorHelixFunctor is a functor for tracking the particle to the spherical mirror
 */
class MirrorHelixFunctor: public HelixFunctor {
 public:
  /**
   * Constructor that sets up the mirror position and curvature
   * @param MirrorPosition Position of mirror in global coordinates
   * @param Curvature Mirror curvature
   */
  MirrorHelixFunctor(const Vector &MirrorPosition, double Curvature);
  /**
   * We do not need copy constructor
   */
  MirrorHelixFunctor(const MirrorHelixFunctor &Functor) = delete;
  /**
   * Virtual () operator
   * @param Vect Vector position
   */
  virtual double operator()(const Vector &Vect) const override;
 private:
  /**
   * The mirror position
   */
  const Vector m_MirrorPosition;
  /**
   * The mirror radius of curvature
   */
  const double m_Curvature;
};

#endif
