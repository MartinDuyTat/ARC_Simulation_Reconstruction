// Martin Duy Tat 29th April 2022

#include<string>
#include<iostream>
#include"TGraph.h"
#include"TCanvas.h"
#include"TLegend.h"
#include"SiPM.h"
#include"Photon.h"

SiPM::SiPM(): m_DetectorSizeX(0.05),
	      m_DetectorSizeY(0.05),
	      m_DetectorPositionX(0.0),
	      m_DetectorPositionY(0.0) {
}

void SiPM::AddPhotonHit(const Photon &photon) {
  m_PhotonHits.emplace_back(photon.m_Position.X(), photon.m_Position.Y(), &photon);
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
  TCanvas c("c", "", 800, 800);
  Graph_Aerogel.SetMarkerStyle(kFullDotLarge);
  Graph_Gas.SetMarkerStyle(kFullDotLarge);
  Graph_Aerogel.SetMarkerColor(kBlue);
  Graph_Gas.SetMarkerColor(kRed);
  Graph_Aerogel.SetTitle("Photon hits;x (cm); y (cm)");
  Graph_Aerogel.Draw("AP");
  Graph_Gas.Draw("P SAME");
  TLegend Legend(0.7, 0.8, 0.9, 0.9);
  Legend.AddEntry(&Graph_Aerogel, "Aerogel 5 GeV", "P");
  Legend.AddEntry(&Graph_Gas, "Gas 5 GeV", "P");
  Legend.Draw();
  c.Draw();
  c.SaveAs(Filename.c_str());
}


const std::vector<PhotonHit>& SiPM::GetPhotonHits() const {
  return m_PhotonHits;
}