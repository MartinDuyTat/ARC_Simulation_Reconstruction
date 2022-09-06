// Martin Duy Tat 2nd September 2022

#include"RadiatorArray.h"
#include"Settings.h"

RadiatorArray::RadiatorArray():
  m_NumberMainRowCells(Settings::GetInt("ARCGeometry/CellsPerRow")),
  m_xHexDist(Settings::GetDouble("ARCGeometry/Length")/(2*m_NumberMainRowCells - 1)),
  m_yHexDist(m_xHexDist*TMath::Sqrt(3)/2.0) {
}
