// Martin Duy Tat 29th April 2022

#include<string>
#include<iostream>
#include<memory>
#include"TGraph.h"
#include"TCanvas.h"
#include"TLegend.h"
#include"TRandom.h"
#include"TBox.h"
#include"TAxis.h"
#include"TPad.h"
#include"Math/Interpolator.h"
#include"SiPM.h"
#include"Photon.h"
#include"Settings.h"

SiPM::SiPM(double xPosition, double yPosition): m_DetectorSizeX(0.08),
						m_DetectorSizeY(0.08),
						m_DetectorPositionX(xPosition),
						m_DetectorPositionY(yPosition),
						m_MaxPDE(0.432),
						m_Interpolator(std::make_unique<Interpolator>(16, InterpolationType::kAKIMA)) {
  m_Interpolator->SetData(16, GetPDEWavelengths().data(), GetMeasuredPDE().data());
}

void SiPM::AddPhotonHit(Photon &photon) {
  m_PhotonHits.emplace_back(photon.m_Position.X(), photon.m_Position.Y(), &photon);
  if(photon.m_Status == Photon::Status::AerogelScattered) {
    return;
  }
  const double PhotonLambda = 1239.8/photon.m_Energy;
  if(PhotonLambda < GetPDEWavelengths()[0] || PhotonLambda > GetPDEWavelengths().back()) {
    std::cout << "Warning! Generated photon outside energy range of SiPM!\n";
  }
  if(!IsDetectorHit(photon)) {
    photon.m_Status = Photon::Status::DetectorMiss;
    return;
  }
  const double RandomNumber = gRandom->Uniform(0.0, m_MaxPDE);
  const double Efficiency = m_Interpolator->Eval(PhotonLambda)*0.90*0.80;
  if(RandomNumber <= Efficiency) {
    photon.m_Status = Photon::Status::DetectorHit;
  } else {
    photon.m_Status = Photon::Status::EfficiencyMiss;
  }
}

void SiPM::PlotHits(const std::string &Filename) const {
  std::vector<double> xHits_Aerogel, yHits_Aerogel, xHits_Gas, yHits_Gas;
  for(const auto &PhotonHit : m_PhotonHits) {
    if(PhotonHit.m_Photon->m_Radiator == Photon::Radiator::Aerogel) {
      xHits_Aerogel.push_back(PhotonHit.x*100.0);
      yHits_Aerogel.push_back(PhotonHit.y*100.0);
    } else if(PhotonHit.m_Photon->m_Radiator == Photon::Radiator::Gas) {
      xHits_Gas.push_back(PhotonHit.x*100.0);
      yHits_Gas.push_back(PhotonHit.y*100.0);
    } else {
      std::cout << "Warning! Detector hit with unknown origin\n";
    }
  }
  TGraph Graph_Aerogel(xHits_Aerogel.size(), xHits_Aerogel.data(), yHits_Aerogel.data());
  TGraph Graph_Gas(xHits_Gas.size(), xHits_Gas.data(), yHits_Gas.data());
  TCanvas c("c", "", 1200, 1200);
  Graph_Aerogel.SetMarkerStyle(kFullDotLarge);
  Graph_Gas.SetMarkerStyle(kFullDotLarge);
  Graph_Aerogel.SetMarkerColor(kBlue);
  Graph_Gas.SetMarkerColor(kRed);
  Graph_Aerogel.SetTitle("Photon hits;x (cm); y (cm)");
  Graph_Aerogel.Draw("AP");
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  Graph_Aerogel.GetXaxis()->SetLimits((m_DetectorPositionX - m_DetectorSizeX/2.0)*100.0,
				      (m_DetectorPositionX + m_DetectorSizeX/2.0)*100.0);
  Graph_Aerogel.SetMinimum((m_DetectorPositionY - m_DetectorSizeY/2.0)*100.0);
  Graph_Aerogel.SetMaximum((m_DetectorPositionY + m_DetectorSizeY/2.0)*100.0);
  Graph_Aerogel.Draw("AP");
  Graph_Gas.Draw("P SAME");
  TLegend Legend(0.7, 0.8, 0.9, 0.9);
  const int Momentum = Settings::GetInt("Particle/Momentum");
  Legend.AddEntry(&Graph_Aerogel, ("Aerogel " + std::to_string(Momentum) + " GeV").c_str(), "P");
  Legend.AddEntry(&Graph_Gas, ("Gas " + std::to_string(Momentum) + " GeV").c_str(), "P");
  Legend.Draw();
  c.Draw();
  c.SaveAs(Filename.c_str());
}

std::unique_ptr<TObject> SiPM::DrawSiPM(const Vector &RadiatorPosition) const {
  const double RadiatorZ = RadiatorPosition.Z();
  const double DetectorPositionX = m_DetectorPositionX + RadiatorPosition.X();
  TBox Box(DetectorPositionX - m_DetectorSizeX/2.0, RadiatorZ - 0.002,
	   DetectorPositionX + m_DetectorSizeX/2.0, RadiatorZ - 0.001);
  Box.SetFillColor(kBlack);
  return std::make_unique<TBox>(Box);
}

bool SiPM::IsDetectorHit(const Photon &photon) const {
  if(TMath::Abs(photon.m_Position.X() - m_DetectorPositionX) > m_DetectorSizeX/2.0) {
    return false;
  } else if(TMath::Abs(photon.m_Position.Y() - m_DetectorPositionY) > m_DetectorSizeY/2.0) {
    return false;
  } else {
    return true;
  }
}

const std::vector<PhotonHit>& SiPM::GetPhotonHits() const {
  return m_PhotonHits;
}

constexpr std::array<double, 16> SiPM::GetPDEWavelengths() const {
  constexpr std::array<double, 16> Lambda{287.0, 299.0, 322.0, 340.0, 367.0, 392.0, 402.0, 413.0,
                                          422.0, 437.0, 452.0, 467.0, 503.0, 590.0, 700.0, 800.0};
  return Lambda;
}

constexpr std::array<double, 16> SiPM::GetMeasuredPDE() const {
  constexpr std::array<double, 16> PDE{0.20, 0.41, 0.46, 0.47, 0.52, 0.59, 0.60, 0.60,
                                       0.59, 0.57, 0.56, 0.53, 0.51, 0.40, 0.26, 0.13};
  return PDE;
}
