/* AUTORIGHTS
Copyright (C) 2006-2018  Ilker R. Capoglu and Di Zhang

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef READ_POINTSOURCES_H
#define READ_POINTSOURCES_H

//use the libconfig library
#include "config/cfg_utils.h"

//only the declaration of Cpointsources needed: use forward declaration
class Cpointsources;
//only the declaration of Cwfs needed: use forward declaration
class Cwfs;

void read_pointsources(Cpointsources &PointSources, const Config& fdtdconfig, const Config& validsettings, const Cwfs& Waveforms);		//modifies the PointSources object

#endif
