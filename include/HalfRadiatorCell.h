// Martin Duy Tat 27th August 2022
/**
 * HalfRadiatorCell inherits from RadiatorCell and is more or less identical, but the cell is only the left half
 */

#ifndef HALFRADIATORCELL
#define HAlFRADIATORCELL

#include"BarrelRadiatorCell.h"

class HalfRadiatorCell: public BarrelRadiatorCell {
 public:
  /**
   * Default constructor that simply calls the parent constructor
   */
  HalfRadiatorCell(std::size_t CellColumnNumber, std::size_t CellRowNumber, double HexagonSize);
  /**
   * Function that checks if the position is inside the radiator
   */
  virtual bool IsInsideCell(const Vector &Position) const override;
  /**
   * Draw radiator geometry
   */
  virtual std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
  DrawRadiatorGeometry() const override;
};

#endif
