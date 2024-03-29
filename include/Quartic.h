/***************************************************************************
 *   Copyright (C) 2016 by Саша Миленковић                                 *
 *   sasa.milenkovic.xyz@gmail.com                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *   ( http://www.gnu.org/licenses/gpl-3.0.en.html )                       *
 *   *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef QUARTIC_H_INCLUDED
#define QUARTIC_H_INCLUDED

#include<complex>
#include<array>

namespace Quartic {

  const double eps=1e-12;

  //---------------------------------------------------------------------------
  // x - array of size 3
  // In case 3 real roots: => x[0], x[1], x[2], return 3
  //         2 real roots: x[0], x[1],          return 2
  //         1 real root : x[0], x[1] ± i*x[2], return 1
  unsigned int solveP3(double* x, double a, double b, double c);

  //---------------------------------------------------------------------------
  // Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
  std::array<std::complex<double>, 4> solve_quartic(double a, double b, double c, double d);

}

#endif // QUARTIC_H_INCLUDED
