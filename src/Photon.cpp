// Martin Duy Tat 29th April 2022

#include<memory>
#include<vector>
#include<utility>
#include<string>
#include"TLine.h"
#include"Photon.h"
#include"RadiatorCell.h"
#include"BarrelRadiatorCell.h"
#include"Utilities.h"

Photon::Photon(const Vector &Position,
	       const Vector &Direction,
	       double Energy,
	       double CosCherenkovAngle,
	       Radiator radiator,
	       const RadiatorCell *radiatorCell):
               m_Position(Position),
	       m_EmissionPoint(Position),
	       m_Direction(Direction),
	       m_Energy(Energy),
	       m_Radiator(radiator),
	       m_CosCherenkovAngle(CosCherenkovAngle),
	       m_Status(Status::Emitted),
	       m_AerogelTravelDistance(0.0),
	       m_MirrorHitPosition(nullptr),
               m_RadiatorCell(radiatorCell) {
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
Photon::DrawPhotonPath() const {
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> PhotonLine;
  const auto RadiatorPosition = m_RadiatorCell->GetRadiatorPosition();
  const auto EmissionPoint = Utilities::SwapXZForEndCap(m_RadiatorCell,
							m_EmissionPoint
							+ RadiatorPosition);
  const auto Position = Utilities::SwapXZForEndCap(m_RadiatorCell,
						   m_Position + RadiatorPosition);
  if(!m_MirrorHitPosition) {
    TLine PhotonLine1(EmissionPoint.X(),
		      EmissionPoint.Z(),
		      Position.X(),
		      Position.Z());
    if(m_Status == Status::MirrorMiss) {
      PhotonLine1.SetLineColor(7);
    } else {
      PhotonLine1.SetLineColor(kGreen);
    }
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine1), ""));
  } else {
    const auto MirrorHitPosition = Utilities::SwapXZForEndCap(m_RadiatorCell,
							      *m_MirrorHitPosition 
							      + RadiatorPosition);
    TLine PhotonLine1(EmissionPoint.X(),
		      EmissionPoint.Z(),
		      MirrorHitPosition.X(),
		      MirrorHitPosition.Z());
    TLine PhotonLine2(MirrorHitPosition.X(),
		      MirrorHitPosition.Z(),
		      Position.X(),
		      Position.Z());
    if(m_Status != Status::DetectorHit) {
      PhotonLine1.SetLineColor(6);
      PhotonLine2.SetLineColor(6);
    } else {
      PhotonLine1.SetLineColor(kBlue);
      PhotonLine2.SetLineColor(kBlue);
    }
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine1), ""));
    PhotonLine.push_back(std::make_pair(std::make_unique<TLine>(PhotonLine2), ""));
  }
  return PhotonLine;
}

Photon::Photon(const Photon &photon):
  m_Position(photon.m_Position),
  m_EmissionPoint(photon.m_EmissionPoint),
  m_Direction(photon.m_Direction),
  m_Energy(photon.m_Energy),
  m_Radiator(photon.m_Radiator),
  m_CosCherenkovAngle(photon.m_CosCherenkovAngle),
  m_Status(photon.m_Status),
  m_AerogelTravelDistance(photon.m_AerogelTravelDistance),
  m_RadiatorCell(photon.m_RadiatorCell) {
  if(photon.m_MirrorHitPosition) {
    m_MirrorHitPosition = std::make_unique<Vector>(*photon.m_MirrorHitPosition);
  } else {
    m_MirrorHitPosition = nullptr;
  }
}
