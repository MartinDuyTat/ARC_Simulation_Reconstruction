// Martin Duy Tat 29th April 2022

#include<memory>
#include<vector>
#include<utility>
#include<string>
#include"TLine.h"
#include"Photon.h"
#include"RadiatorCell.h"

Photon::Photon(const Vector &Position,
	       const Vector &Direction,
	       double Energy,
	       double CherenkovAngle): m_Position(Position),
				       m_EmissionPoint(Position),
				       m_Direction(Direction),
				       m_Energy(Energy),
				       m_Radiator(Radiator::Unknown),
				       m_CherenkovAngle(CherenkovAngle),
				       m_Status(Status::Emitted),
				       m_MirrorHitPosition(nullptr) {
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>> Photon::DrawPhotonPath(const RadiatorCell &radiatorCell) const {
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> PhotonLine;
  const auto EmissionPoint = m_EmissionPoint + radiatorCell.GetRadiatorPosition();
  const auto Position = m_Position + radiatorCell.GetRadiatorPosition();
  if(!m_MirrorHitPosition) {
    TLine PhotonLine1(EmissionPoint.X(), EmissionPoint.Z(), Position.X(), Position.Z());
    if(m_Status == Status::MissedTheta) {
      PhotonLine1.SetLineColor(7);
    } else if(m_Status == Status::MissedPhi) {
      PhotonLine1.SetLineColor(5);
    } else {
      PhotonLine1.SetLineColor(kGreen);
    }
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine1), ""));
  } else {
    const auto MirrorHitPosition = *m_MirrorHitPosition + radiatorCell.GetRadiatorPosition();
    TLine PhotonLine1(EmissionPoint.X(), EmissionPoint.Z(), MirrorHitPosition.X(), MirrorHitPosition.Z());
    TLine PhotonLine2(MirrorHitPosition.X(), MirrorHitPosition.Z(), Position.X(), Position.Z());
    PhotonLine1.SetLineColor(kBlue);
    PhotonLine2.SetLineColor(kBlue);
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine1), ""));
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine2), ""));
  }
  return PhotonLine;
}
