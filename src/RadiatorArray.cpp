// Martin Duy Tat 23rd June 2022

#include<vector>
#include<stdexcept>
#include"Math/Vector3Dfwd.h"
#include"RadiatorArray.h"
#include"RadiatorCell.h"
#include"Settings.h"

RadiatorArray::RadiatorArray(): m_FullArray(Settings::GetBool("General/FullArray")),
                                m_NumberCells(Settings::GetInt("ARCGeometry/ThetaCells")) {
  if(m_FullArray) {
    if(m_NumberCells % 2 == 1) {
      throw std::invalid_argument("Number of cells in theta direction must be even");
    }
    m_Cells.reserve(m_NumberCells);
    for(int i = -m_NumberCells/2; i <= m_NumberCells/2; i++) {
      if(i == 0) {
	continue;
      }
      m_Cells.emplace_back(i);
    }
  } else {
    m_Cells.emplace_back(0);
  }
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>> RadiatorArray::DrawRadiatorArray() const {
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> RadiatorArrayObjects;
  for(const auto &Cell : m_Cells) {
    auto RadiatorObjects = Cell.DrawRadiatorGeometry();
    RadiatorArrayObjects.insert(RadiatorArrayObjects.end(),
				std::make_move_iterator(RadiatorObjects.begin()),
				std::make_move_iterator(RadiatorObjects.end()));
  }
  return RadiatorArrayObjects;
}

RadiatorIter RadiatorArray::operator[](int i) {
  if(i == 0 && m_FullArray) {
    throw std::invalid_argument("Radiator 0 does not exist");
  }
  if(i != 0 && !m_FullArray) {
    throw std::invalid_argument("Only a single radiator cell is present, cannot have radiator cell index " + i);
  }
  if(m_FullArray) {
    if(i < 0) {
      return m_Cells.begin() + (i + m_NumberCells/2);
    } else {
      return m_Cells.begin() + (i - 1 + m_NumberCells/2);
    }
  } else {
    return m_Cells.begin();
  }
}

RadiatorIter RadiatorArray::WhichRadiator(const Vector &Position) {
  const double ThetaLength = m_Cells[0].GetThetaLength();
  const double AbsPosition = TMath::Abs(Position.X());
  if(m_FullArray) {
    const int Sign = Position.X() >= 0.0 ? +1 : -1;
    return (*this)[Sign*(static_cast<int>(AbsPosition/ThetaLength) + 1)];
  } else {
    if(AbsPosition > ThetaLength/2.0) {
      throw std::runtime_error("Particle outside of theta range");
    } else {
      return m_Cells.begin();
    }
  }
}
