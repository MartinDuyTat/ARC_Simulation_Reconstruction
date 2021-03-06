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
   */
  enum class Status{Emitted, MirrorHit, DetectorHit, DetectorMiss, MissedTheta, MissedPhi, MissedDetectorPlane};
  /**
   * Construct a photon with position, direction and energy
   * @param Position Position vector
   * @param Direction Direction vector
   * @param Energy Photon energy
   * @param CherenkovAngle Cherenkov angle
   * @param radiator Aerogel or gas radiator where this photon was emitted
   * @param radiatorCell The radiator cell where this photon was emitted
   */
  Photon(const Vector &Position,
	 const Vector &Direction,
	 double Energy,
	 double CherenkovAngle,
	 Radiator radiator,
	 RadiatorCell *radiatorCell);
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
  const double m_CherenkovAngle;
  /**
   * Flag specifying the status of the photon
   */
  Status m_Status;
  /**
   * Hit position of the mirror, or nullptr if the photon missed the mirror
   */
  std::unique_ptr<Vector> m_MirrorHitPosition;
  /**
   * Pointer to the radiator cell that this photon is inside
   */
  RadiatorCell *m_RadiatorCell;
};

#endif
