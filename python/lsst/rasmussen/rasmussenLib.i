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
%include "file.i"

%lsst_exceptions()

%{
#define PY_ARRAY_UNIQUE_SYMBOL LSST_RASMUSSEN_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"

#include "lsst/pex/logging.h"
#include "lsst/afw.h"
#include "lsst/meas/algorithms.h"
%}

%import "lsst/meas/algorithms/algorithmsLib.i"

%include "ndarray.i"
%init %{
    import_array();
%}

%declareNumPyConverters(ndarray::Array<int,2,2>);

%shared_ptr(lsst::rasmussen::Fe55Control)
%shared_ptr(data_str)
%shared_ptr(lsst::rasmussen::Event)

%{
#include "lsst/rasmussen/Event.h"
#include "lsst/rasmussen/fe55.h"
#include "lsst/rasmussen/tables.h"
%}

%include "lsst/rasmussen/rv.h"
%include "lsst/rasmussen/Event.h"
%include "lsst/rasmussen/fe55.h"
%include "lsst/rasmussen/tables.h"

%template(vectorEvent) std::vector<boost::shared_ptr<lsst::rasmussen::Event> >;

%extend lsst::rasmussen::Event {
    %pythoncode {
    def getData(self, *args):
        """Access the int data[9] array (indexing on Event also works: ev[3])"""
        i = args[0]
        if i < 0 or i >= 9:
            raise IndexError("Index %d is out of range 0..8" % i)

        import ctypes
        __data = 9*ctypes.c_float
        __data = __data.from_address(int(self.data))
        return __data[i]

    __getitem__ = getData
    }
}
