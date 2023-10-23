// Martin Duy Tat 16th October 2022
/**
 * ARCVector is a vector that can optionally hold a local coordinate system
 * The local coordinate system is displaced and rotated from the global one
 * Vector operations are performed on the global coordiate system if no optional
 * coordinate system is held, otherwise vector operations are relative to the
 * local coodinate system
 */

#ifndef ARCVECTOR
#define ARCVECTOR

#include<utility>
#include<optional>
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Math/Rotation3D.h"

using Vector = ROOT::Math::XYZVector;
using Rotation3D = ROOT::Math::Rotation3D;

class ARCVector {
 public:
  /**
   * Constructor that initialises the vector at the origin in global coordinates
   */
  ARCVector() = default;
  /**
   * Constructor that initialises the vector in the global coordinate system
   * @param Vec The global vector
   */
  explicit ARCVector(const Vector &Vec);
  /**
   * Constructor that initialises the vector directly in the local coordinates
   * @param Vec The local vector
   * @param RadiatorPosition The position of the origin of the local coordinates
   * @param RadiatorRotation The rotation from global to local coordinates
   */
  ARCVector(const Vector &Vec,
	    const Vector &RadiatorPosition,
	    const Rotation3D &RadiatorRotation);
  /**
   * Copy constructor is default
   */
  ARCVector(const ARCVector &Vec) = default;
  /**
   * Assign a local coordinate system
   * @param Position The position of the local coordinate system
   * @param Rotation The rotation from the global to local coordinate system
   */
  void AssignLocalCoordinates(const Vector &Position, const Rotation3D &Rotation);
  /**
   * Returns true if the local coordinate system is defined
   */
  bool IsLocal() const;
  /**
   * Delete local coordinates and turn this vector global again
   */
  void ConvertToGlobal();
  /**
   * Get the local vector
   */
  const Vector& LocalVector() const;
  /**
   * Get the global vector
   */
  Vector GlobalVector() const;
  /**
   * Operator for adding a vector
   */
  ARCVector& operator +=(const Vector &Vec);
  /**
   * Rotate vector in the azimuthal direction
   */
  void MapPhi(double DeltaPhi);
  /**
   * Reflect vector along z direction
   */
  void ReflectZ();
  /**
   * Reflect vector along y direction
   */
  void ReflectY();
  /**
   * Set the x-component, only use this to set the mirror position!
   */
  void SetX(double x);
  /**
   * Set the z-component, only use this to set the mirror position!
   */
  void SetZ(double z);
  /**
   * Set global vector
   */
  void SetGlobalVector(const Vector &Vect);
 private:
  /**
   * The actual vector
   * If no local coordinate system is defined, it's a global vector
   * If a local coordinate system is defined, it's a local vector
   */
  Vector m_Vector;
  /**
   * The optional position and rotation of the local coordinate system
   */
  std::optional<std::pair<Vector, Rotation3D>> m_CoordinateTransform;
  /**
   * Helper function that transforms a vector from global to local coordinates
   * @param Vec Vector to be transformed
   */
  void GlobalToLocal(Vector &Vec) const;
  /**
   * Helper function that transforms a vector from local to global coordinates
   * @param Vec Vector to be transformed
   */
  void LocalToGlobal(Vector &Vec) const;
  
};

#endif
