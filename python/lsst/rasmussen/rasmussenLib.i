// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
%define rasmussenLib_DOCSTRING
"
Various swigged-up C++ classes for rasmussen
"
%enddef

%feature("autodoc", "1");
%module(package="rasmussenLib", docstring=rasmussenLib_DOCSTRING) rasmussenLib

%pythonnondynamic;
%naturalvar;  // use const reference typemaps

%include "lsst/p_lsstSwig.i"

%lsst_exceptions()

%{
#include "lsst/pex/logging.h"
#include "lsst/afw.h"
#include "lsst/meas/algorithms.h"
%}

%import "lsst/meas/algorithms/algorithmsLib.i"

%shared_ptr(lsst::rasmussen::Fe55Control)

%{
#include "lsst/rasmussen/Event.h"
#include "lsst/rasmussen/fe55.h"
#include "lsst/rasmussen/tables.h"
%}

%include "lsst/rasmussen/rv.h"
%include "lsst/rasmussen/Event.h"
%include "lsst/rasmussen/fe55.h"
%include "lsst/rasmussen/tables.h"
