// Martin Duy Tat 29th April 2022

#include<memory>
#include"TLine.h"
#include"Photon.h"

Photon::Photon(const Vector &Position,
	       const Vector &Direction,
	       double Energy,
	       double CherenkovAngle): m_Position(Position),
				       m_EmissionPoint(Position),
				       m_Direction(Direction),
				       m_Energy(Energy),
				       m_Radiator(Radiator::Unknown),
				       m_CherenkovAngle(CherenkovAngle),
				       m_MirrorHit(false) {
}

std::unique_ptr<TLine> Photon::DrawPhotonPath() const {
  TLine PhotonLine(m_EmissionPoint.X(), m_EmissionPoint.Z(),
		   m_Position.X(), m_Position.Z());
  PhotonLine.SetLineColor(kBlue);
  return std::make_unique<TLine>(PhotonLine);
}
