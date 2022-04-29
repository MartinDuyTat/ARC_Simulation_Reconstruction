// Martin Duy Tat 29th April 2022
/**
 * A Photon has a position, direction and energy
 */

#ifndef PHOTON
#define PHOTON

#include"Vector3Dfwd.h"

using Vector = ROOT::Math::XYZVector;

struct Photon {
  /**
   * Construct a photon with position, direction and energy
   * @param Position Position vector
   * @param Direction Direction vector
   * @param Energy Photon energy
   */
  Photon(const Vector &Position, const Vector &Direction, double Energy);
  /**
   * Photon position in local cell coordinates, in m
   */
  Vector m_Position;
  /**
   * Photon direction vector
   */
  Vector m_Direction;
  /**
   * Photon energy
   */
  double m_Energy;
};

#endif
