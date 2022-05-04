// Martin Duy Tat 3rd May 2022

#include"ReconstructedPhoton.h"

ReconstructedPhoton::ReconstructedPhoton(const Photon &photon): m_Photon(&photon),
								m_CherenkovAngle_TrueEmissionPoint(0.0),
								m_CherenkovAngle(0.0) {
}
