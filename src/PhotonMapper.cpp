// Martin Duy Tat 1st May 2022

#include<algorithm>
#include"TMath.h"
#include"PhotonMapper.h"
#include"RadiatorCell.h"
#include"Photon.h"

namespace PhotonMapper {

  double PhotonMirrorDistance(const Photon &photon, const RadiatorCell &radiatorCell) {
    auto MirrorCentre = radiatorCell.GetMirrorCentre();
    const double R = radiatorCell.GetMirrorCurvature();
    const double a = 1.0;
    const double b = -2*photon.m_Direction.Dot(MirrorCentre - photon.m_Position);
    const double c = (MirrorCentre - photon.m_Position).Mag2() - R*R;
    const double s1 = (-b + TMath::Sqrt(b*b - 4.0*a*c))/(2.0*a);
    const double s2 = (-b - TMath::Sqrt(b*b - 4.0*a*c))/(2.0*a);
    return std::max(s1, s2);
  }

  void TracePhotonToMirror(Photon &photon, const RadiatorCell &radiatorCell) {
    const double s = PhotonMirrorDistance(photon, radiatorCell);
    const Vector Mirror = photon.m_Position + s*photon.m_Direction;
    const Vector Normal = (Mirror - radiatorCell.GetMirrorCentre())/radiatorCell.GetMirrorCurvature();
    const Vector Vout = photon.m_Direction - 2*photon.m_Direction.Dot(Normal)*Normal;
    photon.m_Position = Mirror;
    photon.m_Direction = Vout;
    // TODO: Check that photon hits mirror
    if(true) {
      photon.m_MirrorHit = true;
    }
  }

  void TracePhotonToDetector(Photon &photon, RadiatorCell &radiatorCell) {
    photon.m_Position += (photon.m_Position.Dot(Vector(0.0, 0.0, 1.0))/TMath::Abs(photon.m_Direction.Z()))*photon.m_Direction;
    radiatorCell.m_Detector.AddPhotonHit(photon);
  }

  void TracePhoton(Photon &photon, RadiatorCell &radiatorCell) {
    TracePhotonToMirror(photon, radiatorCell);
    TracePhotonToDetector(photon, radiatorCell);
  }

}
