// Martin Duy Tat 29th September 2022
/**
 * Particle is an abstract class that stores the position and radiator cell that the particle is in
 */

#ifndef PARTICLE
#define PARTICLE

#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"

using Vector = ROOT::Math::XYZVector;

class RadiatorCell;

class Particle {
 public:
  /**
   * Enum with the two coordinate systems used
   * GlobalDetector: z along symmetry axis, x is up, y is out of the plane (only used to generate charged tracks)
   * LocalRadiator: z axis in the radial direction, x in the theta direction and y in the phi direction
   */
  enum class CoordinateSystem{GlobalDetector, LocalRadiator};
  /**
   * Constructor that stores the position and radiator cell
   * @param Position The particle position
   * @param coordinateSystem The coordinate system that this particle is in
   * @param radiatorCell The radiator cell the particle belong to
   */
  Particle(const Vector &Position,
	   CoordinateSystem coordinateSystem,
	   const RadiatorCell *radiatorCell = nullptr);
  /**
   * Need virtual constructor since it's a virtual class
   */
  ~Particle() = default;
  /**
   * Convert to local radiator coordinates
   */
  virtual void ConvertToRadiatorCoordinates();
  /**
   * Convert back to global coordinates
   */
  virtual void ConvertBackToGlobalCoordinates();
  /**
   * Get particle position
   */
  const Vector& GetPosition() const;
  /**
   * Rotate around the phi direction to map particle to a valid radiator cell
   * @param DeltaPhi Angle that is rotated in phi
   */
  virtual void MapPhi(double DeltaPhi);
  /**
   * Reflect everything in the z-direction since radiator cells are symmetric
   */
  virtual void ReflectZ();
  /**
   * Reflect everything in the y-direction
   */
  virtual void ReflectY();
  /**
   * Set the rotation angle in phi relative to the global detector coordinates
   */
  void SetPhiRotated(double Phi);
  /**
   * Get a pointer to the radiator cell
   */
  const RadiatorCell* GetRadiatorCell() const;
  /**
   * Check if particle is at the radiator so that we can search for the radiator cell
   */
  virtual bool IsAtRadiator() const = 0;
 protected:
  /**
   * Particle position, in m
   */
  Vector m_Position;
  /**
   * Pointer to the radiator cell that track enters
   */
  const RadiatorCell *m_RadiatorCell;
  /**
   * Helper function to swap x and z directions (so that z points up in the cell)
   */
  void SwapXZ(Vector &Vec) const;
  /**
   * Helper function to swap x and z directions of all vectors (where relevant)
   */
  virtual void SwapXZ();
  /**
   * Helper function to change coordinate origin
   */
  virtual void ChangeCoordinateOrigin(const Vector &Shift);
  /**
   * Which coordinate system the momentum and position is given in
   */
  CoordinateSystem m_CoordinateSystem;
  /**
   * The rotation in phi relative to the global detector coordinates
   */
  double m_PhiRotated;
};

#endif
