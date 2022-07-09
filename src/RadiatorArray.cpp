// Martin Duy Tat 23rd June 2022

#include<vector>
#include<stdexcept>
#include<string>
#include"Math/Vector3Dfwd.h"
#include"RadiatorArray.h"
#include"RadiatorCell.h"
#include"Settings.h"

RadiatorArray::RadiatorArray(): m_FullArray(Settings::GetBool("General/FullArray")),
                                m_NumberMainRowCells(Settings::GetInt("ARCGeometry/CellsPerRow")),
				m_xHexDist(Settings::GetDouble("ARCGeometry/Length")/(2*m_NumberMainRowCells + 1)),
				m_yHexDist(m_xHexDist/TMath::Sqrt(3)) {
  if(m_FullArray) {
    m_Cells.reserve(2*m_NumberMainRowCells);
    for(int i = 0; i < m_NumberMainRowCells; i++) {
      // Main row
      m_Cells.emplace_back(1, i, m_xHexDist);
    }
    for(int i = 0; i < m_NumberMainRowCells; i++) {
      // Upper row
      m_Cells.emplace_back(2, i, m_xHexDist);
    }
  } else {
    m_Cells.emplace_back(0, 0, m_xHexDist);
  }
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
RadiatorArray::DrawRadiatorArray() const {
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> RadiatorArrayObjects;
  for(const auto &Cell : m_Cells) {
    auto RadiatorObjects = Cell.DrawRadiatorGeometry();
    RadiatorArrayObjects.insert(RadiatorArrayObjects.end(),
				std::make_move_iterator(RadiatorObjects.begin()),
				std::make_move_iterator(RadiatorObjects.end()));
  }
  return RadiatorArrayObjects;
}

RadiatorIter RadiatorArray::operator()(int i, int j) {
  if(m_FullArray) {
    if(j == 0 || (i != 1 && i != 2)) {
      throw std::invalid_argument("Radiator (" + std::to_string(i) + ", " + std::to_string(j) + ") does not exist");
    } else {
      if(i == 1) {
	return m_Cells.begin() + j;
      } else {
	return m_Cells.begin() + m_NumberMainRowCells + j;
      }
    }
  } else {
    if(i != 0 && j != 0) {
      throw std::invalid_argument("Cannot have radiator cell index (" + std::to_string(i) + ", " + std::to_string(j) + " )");
    } else {
      return m_Cells.begin();
    }
  }
}

RadiatorIter RadiatorArray::WhichRadiator(const Vector &Position) {
  if(m_FullArray) {
    const int yIndex = std::round(Position.Y()/m_yHexDist) + 1;
    const bool MainRow = yIndex % 2 == 0;
    const double xDist = MainRow ? Position.X() : (Position.X() + m_xHexDist/2.0);
    int xIndex = std::round(xDist/m_xHexDist);
    RadiatorIter ThisRadiator = (*this)(xIndex, yIndex);
    double Distance2 = (Position - ThisRadiator->GetRadiatorPosition()).Mag2();
    std::vector<std::pair<int, int>> Neighbours{
      {xIndex, yIndex + 1},
      {xIndex, yIndex - 1},
      {xIndex + 1, yIndex},
      {xIndex + 1, yIndex + 1},
      {xIndex - 1, yIndex},
      {xIndex - 1, yIndex + 1}};
    for(const auto &Neighbour : Neighbours) {
      RadiatorIter NeighbourRadiator = (*this)(Neighbour.first, Neighbour.second);
      double NewDistance2 = (Position - NeighbourRadiator->GetRadiatorPosition()).Mag2();
      if(NewDistance2 < Distance2) {
	Distance2 = NewDistance2;
	ThisRadiator = NeighbourRadiator;
      }
    }
    return ThisRadiator;
  } else {
    return m_Cells.begin();
  }
}
