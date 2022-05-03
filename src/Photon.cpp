// Martin Duy Tat 29th April 2022

#include"Photon.h"

Photon::Photon(const Vector &Position, const Vector &Direction, double Energy, double CherenkovAngle): m_Position(Position),
												       m_EmissionPoint(Position),
												       m_Direction(Direction),
												       m_Energy(Energy),
												       m_Radiator(Radiator::Unknown),
												       m_CherenkovAngle(CherenkovAngle) {
}
