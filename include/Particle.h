// Martin Duy Tat 29th September 2022
/**
 * Particle is an abstract class that stores the position and radiator cell that the particle is in
 */

#ifndef PARTICLE
#define PARTICLE

#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"ARCVector.h"

using Vector = ROOT::Math::XYZVector;

class RadiatorCell;
class RadiatorArray;

class Particle {
 public:
  /**
   * Constructor that stores the position and radiator cell
   * @param Position The particle position
   * @param radiatorCell The radiator cell the particle belong to
   */
  Particle(const Vector &Position,
	   const RadiatorCell *radiatorCell = nullptr);
  /**
   * Need virtual constructor since it's a virtual class
   */
  ~Particle() = default;
  /**
   * Find the correct radiator that the particle goes through
   */
  bool FindRadiator(const RadiatorArray &radiatorArray);
  /**
   * Convert to local radiator coordinates
   * @param radiatorCell The radiator cell that we want the local coordinates of
   */
  virtual void ConvertToRadiatorCoordinates();
  /**
   * Convert back to global coordinates
   */
  virtual void ConvertBackToGlobalCoordinates();
  /**
   * Get particle position
   */
  const ARCVector& GetPosition() const;
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
   * Get a pointer to the radiator cell
   */
  const RadiatorCell* GetRadiatorCell() const;
  /**
   * Check if particle is at the radiator so that we can search for the radiator cell
   */
  virtual bool IsAtRadiator() const = 0;
  /**
   * Get the column number of the radiator
   */
  std::size_t GetRadiatorColumnNumber() const;
  /**
   * Get the row number of the radiator
   */
  std::size_t GetRadiatorRowNumber() const;
 protected:
  /**
   * Particle position, in m
   */
  ARCVector m_Position;
  /**
   * Pointer to the radiator cell that track enters
   */
  const RadiatorCell *m_RadiatorCell;
};

#endif
