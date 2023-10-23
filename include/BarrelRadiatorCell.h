// Martin Duy Tat 6th September 2022
/**
 * BarrelRadiatorCell inherits from RadiatorCell and contains the position of the radiator
 */

#ifndef BARRELRADIATORCELL
#define BARRELRADIATORCELL

#include<vector>
#include<memory>
#include<string>
#include"TObject.h"
#include"Math/Rotation3D.h"
#include"RadiatorCell.h"
#include"ARCVector.h"

class BarrelRadiatorCell: public RadiatorCell {
 public:
  /**
   * Constructor that sets up the geometry and a position of the barrel cell
   * @param CellColumnNumber The column number of this cell
   * @param CellRowNumber The row number of this cell
   * @param HexagonSize The length between two opposide hexagon edges (not points)
   */
  BarrelRadiatorCell(std::size_t CellColumnNumber,
		     std::size_t CellRowNumber,
		     double HexagonSize);
  /**
   * Delete copy constructor
   */
  BarrelRadiatorCell(const BarrelRadiatorCell &radiatorCell) = delete;
  /**
   * Function that checks if the position is inside the radiator
   */
  virtual bool IsInsideCell(const ARCVector &Position) const override;
  /**
   * Draw radiator geometry
   */
  virtual std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
  DrawRadiatorGeometry() const override;
  /**
   * Get the reverse rotation in phi, for plotting purposes
   */
  virtual RotationZ ReversePhiRotation() const override;
  /**
   * Function that checks if the detector is contained inside the radiator cell
   */
  virtual bool IsDetectorInsideCell() const override;
  /**
   * Checks if the cell is at the edge (near the end cap)
   */
  virtual bool IsEdgeCell() const override;
 protected:
  /**
   * Radius of the barrel where the entrance window is
   */
  const double m_BarrelRadius;
 private:
  /**
   * Helper function to get cell position based on cell number
   */
  Vector GetCellPosition(std::size_t CellColumnNumber,
			 std::size_t CellRowNumber,
			 double HexagonSize);
  /**
   * Helper function to determine the rotation from global to local coordinates
   */
  Rotation3D GetCellOrientation(double HexagonSize,
				std::size_t CellRowNumber);
};

#endif
