// Martin Duy Tat 23rd June 2022
/**
 * RadiatorArray contains the full ARC geometry of radiators
 */

#ifndef RADIATORARRAY
#define RADIATORARRAY

#include<vector>
#include<memory>
#include<utility>
#include<string>
#include"TObject.h"
#include"Math/Vector3Dfwd.h"
#include"RadiatorCell.h"

using Vector = ROOT::Math::XYZVector;

class RadiatorArray {
 public:
  /**
   * Constructor that sets up the geometry
   */
  RadiatorArray();
  /**
   * Draw all the radiator cells
   */
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> DrawRadiatorArray() const;
  /**
   * Operator overload to get the individual radiator cells
   */
  RadiatorCell& operator[](int i);
  /**
   * Check which radiator the particle goes through
   * It is assumed that the input vector is right at the bottom surface of the radiator plane
   */
  RadiatorCell* WhichRadiator(const Vector &Position);
 private:
  /**
   * Vector containing all the radiator cells
   */
  std::vector<std::unique_ptr<RadiatorCell>> m_Cells;
  /**
   * Flag that is true if the full geometry is included, otherwise a single cell in the middle is used
   */
  const bool m_FullArray;
  /**
   * Number of cells in theta direction
   */
  const int m_NumberCells;
};

#endif
