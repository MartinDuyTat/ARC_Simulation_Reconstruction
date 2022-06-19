// Martin Duy Tat 29th April 2022
/**
 * A Photon has a position, direction and energy
 */

#ifndef PHOTON
#define PHOTON

#include<memory>
#include"TLine.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"

using Vector = ROOT::Math::XYZVector;

struct Photon {
  /**
   * Enum class classifying which radiator the photon was emitted from
   */
  enum class Radiator{Aerogel, Gas, Unknown};
  /**
   * Construct a photon with position, direction and energy
   * @param Position Position vector
   * @param Direction Direction vector
   * @param Energy Photon energy
   * @param CherenkovAngle Cherenkov angle
   */
  Photon(const Vector &Position, const Vector &Direction, double Energy, double CherenkovAngle);
  /**
   * Draw photon path
   */
  std::unique_ptr<TLine> DrawPhotonPath() const;
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
   * Flag that is true when photon has hit the mirror
   */
  bool m_MirrorHit;
};

#endif
