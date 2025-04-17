//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2020, David Radice <david.radice@psu.edu>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef nuX_M1_MACRO_H
#define nuX_M1_MACRO_H

// #ifndef WARN_FOR_SRC_FIX
// #define WARN_FOR_SRC_FIX
// #endif

#define nuX_M1_SRC_EXPL 1  // explicit RHS
#define nuX_M1_SRC_IMPL 2  // implicit RHS (default)
#define nuX_M1_SRC_BOOST 3 // boost to fluid frame (approximate!)

#define nuX_M1_NGHOST 2

#ifndef nuX_M1_SRC_METHOD
#define nuX_M1_SRC_METHOD nuX_M1_SRC_IMPL
#endif

#define SQ(X) ((X) * (X))

#endif
