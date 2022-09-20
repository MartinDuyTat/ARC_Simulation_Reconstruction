// Martin Duy Tat 2nd September 2022

#include"RadiatorArray.h"
#include"Settings.h"

RadiatorArray::RadiatorArray():
  m_NumberMainRowCells(Settings::GetSizeT("ARCGeometry/CellsPerRow")),
  m_xHexDist(Settings::GetDouble("ARCGeometry/Length")/
	     static_cast<double>(2*m_NumberMainRowCells - 1)),
  m_yHexDist(m_xHexDist*TMath::Sqrt(3)/2.0) {
}

RadiatorCell* RadiatorArray::GetRadiatorCell(std::size_t i, std::size_t j) {
  int Index = FindRadiatorIndex(i, j);
  if(Index < 0) {
    return nullptr;
  } else {
    return m_Cells[static_cast<std::size_t>(Index)].get();
  }
}

const RadiatorCell* RadiatorArray::operator()(std::size_t i, std::size_t j) const {
  int Index = FindRadiatorIndex(i, j);
  if(Index < 0) {
    return nullptr;
  } else {
    return m_Cells[static_cast<std::size_t>(Index)].get();
  }
}
