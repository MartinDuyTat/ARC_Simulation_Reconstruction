// Martin Duy Tat 29th April 2022
/**
 * A Photon has a position, direction and energy
 */

#ifndef PHOTON
#define PHOTON

#include<memory>
#include<vector>
#include<utility>
#include<string>
#include"TLine.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Particle.h"
#include"ARCVector.h"

using Vector = ROOT::Math::XYZVector;

class RadiatorCell;

class Photon: public Particle {
 public:
  /**
   * Enum class classifying which radiator the photon was emitted from
   */
  enum class Radiator{Aerogel, Gas};
  /**
   * Enum class classifying whether the photon has been tracked or if it missed the mirror
   * Emitted: Photon just got emitted inside radiator
   * MirrorHit: Photon is at the mirror
   * MirrorMiss: Photon missed the mirror and ended up outside radiator
   * EfficiencyMiss: Photon lost in mirror reflectivity, detector dead space or PDE
   * DetectorHit: Photon hit the detector
   * DetectorMiss: Photon missed the detector
   * AerogelScattered: Scattered by the aerogel
   * WallMiss: The photon hit the vessel wall on top, so an area without a mirror
   * Backwards: The photon must travel backwards to hit mirror (impossible!)
   */
  enum class Status {
    Emitted,
    MirrorHit,
    MirrorMiss,
    EfficiencyMiss,
    DetectorHit,
    DetectorMiss,
    AerogelScattered,
    WallMiss,
    Backwards,
    OutsideCell,
    OutsideMirrorRadius
  };
  /**
   * Construct a photon with position, direction and energy
   * Input vectors in constructor are in global coordinates
   * All vectors are stored in the local coordinates internally
   * @param Position Position vector
   * @param AssumedPosition The midpoint of the track in the radiator
   * @param Direction Direction vector
   * @param ParticleDirection Unit vector in the direction of the particle
   * @param AssumedParticleDirection Unit vector in the assumed particle direction
   * @param Energy Photon energy
   * @param CosCherenkovAngle Cosine of Cherenkov angle
   * @param radiator Aerogel or gas radiator where this photon was emitted
   * @param radiatorCell The radiator cell where this photon was emitted
   * @param IsTrackDrawn True if the track that emitted this photon is drawn
   * @param Weight The weight of this photon, inverse of the photon multiplier
   */
  Photon(const Vector &Position,
	 const Vector &AssumedPosition,
	 const Vector &Direction,
	 const Vector &ParticleDirection,
	 const Vector &AssumedParticleDirection,
	 double Energy,
	 double CosCherenkovAngle,
	 Radiator radiator,
	 const RadiatorCell *radiatorCell,
	 bool IsTrackDrawn,
	 double Weight);
  /**
   * Need virtual constructor since it's a virtual class
   */
  ~Photon() = default;
  /**
   * We don't need copy constructor
   */
  Photon(const Photon &photon) = delete;
  /**
   * The move constructor is default
   */
  Photon(Photon &&photon) = default;
  /**
   * Draw photon path
   */
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> DrawPhotonPath() const;
  /**
   * Propagate the photon
   * @param Displacement The displace that the photon is propagated
   */
  void PropagatePhoton(const Vector &Displacement);
  /**
   * Get the photon direction
   */
  const Vector& GetDirection() const;
  /**
   * Change the photon direction
   * @param Kick The kick given to the photon to change its direction
   */
  void KickPhoton(const Vector &Kick);
  /**
   * Get the photon status
   */
  Status GetStatus() const;
  /**
   * Get the emission point
   */
  const Vector& GetEmissionPoint(bool TrueEmissionPoint) const;
  /**
   * Update photon status
   */
  void UpdatePhotonStatus(Status status);
  /**
   * Get the photon energy
   */
  double GetEnergy() const;
  /**
   * Add travel distance through aerogel
   */
  void AddAerogelTravelDistance(double Distance);
  /**
   * Get travel distance through aerogel
   */
  double GetAerogelTravelDistance() const;
  /**
   * Get mirror hit position, or nullptr if the photon didn't hit the mirror
   */
  const ARCVector* GetMirrorHitPosition() const;
  /**
   * Register mirror hit position
   */
  void RegisterMirrorHitPosition(const Vector &MirrorHitPosition);
  /**
   * Get the radiator that the photon was generated in
   */
  Radiator GetRadiator() const;
  /**
   * Get the cosine of the Cherenkov angle
   */
  double GetCosCherenkovAngle() const;
  /**
   * Convert to local radiator coordinates
   */
  virtual void ConvertToRadiatorCoordinates() override;
  /**
   * Convert back to global coordinates
   */
  virtual void ConvertBackToGlobalCoordinates() override;
  /**
   * Rotate around the phi direction to map particle to a valid radiator cell
   * @param DeltaPhi Angle that is rotated in phi
   */
  virtual void MapPhi(double DeltaPhi) override;
  /**
   * Reflect everything in the z-direction since radiator cells are symmetric
   */
  virtual void ReflectZ() override;
  /**
   * Reflect everything in the y-direction
   */
  virtual void ReflectY() override;
  /**
   * Check if particle is at the radiator so that we can search for the radiator cell
   */
  virtual bool IsAtRadiator() const override;
  /**
   * Set photon position back to the emission point
   */
  void PutPhotonToEmissionPoint();
  /**
   * Call this function is photon has migrated to a different cell
   */
  void PhotonHasMigrated();
  /**
   * Get flag that is true if photon migrates
   */
  bool HasPhotonMigrated() const;
  /**
   * Get the direction of the particle that emitted this photon
   */
  const Vector& GetParticleDirection(bool TrueEmissionPoint) const;
  /**
   * Get the weight of this photon (inverse of the photon multiplier)
   */
  double GetWeight() const;
 private:
  /**
   * The midpoint of the track inside the radiator, where we assume the photon is emitted
   */
  ARCVector m_AssumedEmissionPoint;
  /**
   * Emission point
   */
  ARCVector m_EmissionPoint;
  /**
   * Photon direction vector
   */
  ARCVector m_Direction;
  /**
   * Direction of the particle that emitted this photon
   */
  ARCVector m_ParticleDirection;
  /**
   * The direction of the particle that emitted this photon at the mid point
   */
  ARCVector m_AssumedParticleDirection;
  /**
   * Photon energy
   */
  const double m_Energy;
  /**
   * Flag specifying which radiator the photon was emitted from
   */
  Radiator m_Radiator;
  /**
   * Cherenkov angle
   */
  const double m_CosCherenkovAngle;
  /**
   * Flag specifying the status of the photon
   */
  Status m_Status;
  /**
   * The total distance travelled through aerogel
   */
  double m_AerogelTravelDistance;
  /**
   * Hit position of the mirror, or nullptr if the photon missed the mirror
   */
  std::unique_ptr<ARCVector> m_MirrorHitPosition;
  /**
   * The weight of this photon, usually the inverse of the photon multiplier
   */
  double m_Weight;
  /**
   * Flag that is true if photon migrates to a neighbouring cell
   */
  bool m_HasMigrated;
  /**
   * Flag that indicates if the particle track that emitted this photon is drawn
   */
  bool m_IsTrackDrawn;
};

#endif
