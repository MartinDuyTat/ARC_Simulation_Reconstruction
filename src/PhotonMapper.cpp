// Martin Duy Tat 1st May 2022

#include<algorithm>
#include"TMath.h"
#include"TRandom.h"
#include"PhotonMapper.h"
#include"Photon.h"
#include"Settings.h"

namespace PhotonMapper {

  double PhotonMirrorDistance(const Photon &photon) {
    auto MirrorCentre = photon.m_RadiatorCell->GetMirrorCentre();
    const double R = photon.m_RadiatorCell->GetMirrorCurvature();
    const double a = 1.0;
    const double b = -2*photon.m_Direction.Dot(MirrorCentre - photon.m_Position);
    const double c = (MirrorCentre - photon.m_Position).Mag2() - R*R;
    const double s1 = (-b + TMath::Sqrt(b*b - 4.0*a*c))/(2.0*a);
    const double s2 = (-b - TMath::Sqrt(b*b - 4.0*a*c))/(2.0*a);
    return std::max(s1, s2);
  }

  void TracePhotonToMirror(Photon &photon) {
    const double s = PhotonMirrorDistance(photon);
    const Vector Mirror = photon.m_Position + s*photon.m_Direction;
    const auto Curvature = photon.m_RadiatorCell->GetMirrorCurvature();
    const Vector Normal = (Mirror - photon.m_RadiatorCell->GetMirrorCentre())/Curvature;
    const Vector Vout = photon.m_Direction - 2*photon.m_Direction.Dot(Normal)*Normal;
    photon.m_Position = Mirror;
    photon.m_Direction = Vout;
    if(photon.m_RadiatorCell->IsInsideCell(photon)) {
      photon.m_Status = Photon::Status::MirrorHit;
      photon.m_MirrorHitPosition = std::make_unique<Vector>(photon.m_Position);
    } else {
      photon.m_Status = Photon::Status::MirrorMiss;
    }
  }

  void TracePhotonToDetector(Photon &photon) {
    const double Slope = TMath::Abs(photon.m_Direction.Z());
    // TODO: Account for photons generated in aerogel
    const double AerogelThickness = photon.m_RadiatorCell->GetAerogelThickness();
    constexpr double T0 = 0.9679;
    constexpr double Clarity = 5.087e10;
    const double Lambda = 1239.8/photon.m_Energy;
    const double Transmission = T0*TMath::Exp(-Clarity*(AerogelThickness/Slope)/TMath::Power(Lambda, 4));
    if(gRandom->Uniform(0.0, 1.0) > Transmission) {
      photon.m_Status = Photon::Status::AerogelScattered;
    }
    photon.m_Position += (photon.m_Position.Dot(Vector(0.0, 0.0, 1.0))/Slope)*photon.m_Direction;
    photon.m_RadiatorCell->GetDetector().AddPhotonHit(photon);
  }

  void TracePhoton(Photon &photon) {
    TracePhotonToMirror(photon);
    if(photon.m_MirrorHitPosition) {
      TracePhotonToDetector(photon);
    }
  }

}
