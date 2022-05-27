// Martin Duy Tat 27th May 2022

#include<string>
#include"TCanvas.h"
#include"TPad.h"
#include"EventDisplay.h"
#include"Settings.h"

void EventDisplay::DrawEventDisplay(const std::string &Filename) {
  TCanvas c("c", "", Settings::GetInt("EventDisplay/CanvasWidth"), Settings::GetInt("EventDisplay/CanvasHeight"));
  gPad->DrawFrame(-Settings::GetDouble("ARCGeometry/Length")/2.0, 0.0, Settings::GetDouble("ARCGeometry/Length")/2.0, Settings::GetDouble("ARCGeometry/Radius"));
  for(auto &Object : m_EventObjects) {
    Object->Draw();
  }
  c.Draw();
  c.SaveAs(Filename.c_str());
}

void EventDisplay::AddObject(std::unique_ptr<TObject> Object) {
  m_EventObjects.push_back(std::move(Object));
}

void EventDisplay::AddObject(std::vector<std::unique_ptr<TObject>> Objects) {
  for(auto &Object : Objects) {
    AddObject(std::move(Object));
  }
}
