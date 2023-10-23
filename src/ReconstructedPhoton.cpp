// Martin Duy Tat 3rd May 2022

#include"ReconstructedPhoton.h"

ReconstructedPhoton::ReconstructedPhoton(const Photon &photon):
  m_Photon(&photon),
  m_CosCherenkovAngle_TrueEmissionPoint(1.0),
  m_CosCherenkovAngle(1.0) {
}
