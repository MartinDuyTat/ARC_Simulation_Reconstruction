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

using Vector = ROOT::Math::XYZVector;

class RadiatorCell;

struct Photon {
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
    Backwards
  };
  /**
   * Construct a photon with position, direction and energy
   * @param Position Position vector
   * @param Direction Direction vector
   * @param Energy Photon energy
   * @param CosCherenkovAngle Cosine of Cherenkov angle
   * @param radiator Aerogel or gas radiator where this photon was emitted
   * @param radiatorCell The radiator cell where this photon was emitted
   */
  Photon(const Vector &Position,
	 const Vector &Direction,
	 double Energy,
	 double CosCherenkovAngle,
	 Radiator radiator,
	 const RadiatorCell *radiatorCell);
  /**
   * Copy constructor
   */
  Photon(const Photon &photon);
  /**
   * Draw photon path
   */
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> DrawPhotonPath() const;
  /**
   * Photon position in local cell coordinates, in m
   */
  Vector m_Position;
  /**
   * Emission point
   */
  const Vector m_EmissionPoint;
  /**
   * Photon direction vector
   */
  Vector m_Direction;
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
  std::unique_ptr<Vector> m_MirrorHitPosition;
  /**
   * Pointer to the radiator cell that this photon is inside
   */
  const RadiatorCell *m_RadiatorCell;
};

#endif
