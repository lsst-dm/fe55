#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Process Fe55 data
"""

import os
import os.path
import sys

import lsst.pex.exceptions
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9

import lsst.rasmussen as ras

import numpy

try:
    type(display)
except NameError:
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def processImage(thresh, fileName, grades=range(9), split=None):
    image = afwImage.ImageF(fileName)
    image -= afwMath.makeStatistics(image, afwMath.MEDIAN).getValue()

    if display:
        ds9.mtv(image, title="bkgd subtracted", frame=0)

    fs = afwDetect.FootprintSet(image, afwDetect.Threshold(thresh))

    events = []
    for foot in fs.getFootprints():
        for peak in foot.getPeaks():
            events.append(ras.Event(image, afwGeom.PointI(peak.getIx(), peak.getIy())))

    print "Found %d events" % (len(events))

    if split is None:
        split = int(0.33*thresh)
    filt = sum([1 << g for g in grades])
    table = ras.HistogramTableGflt(filt, thresh, split)

    size = 1.6                          # half-size of box to draw
    ctypes = dict(zip(range(-1, 8), (ds9.WHITE,
                                     ds9.RED, ds9.GREEN, ds9.BLUE,
                                     ds9.CYAN, ds9.MAGENTA, ds9.YELLOW,
                                     ds9.RED, ds9.GREEN)))
    with ds9.Buffering():
        for ev in events:
            if not table.process_event(ev):
                if display:
                    ds9.dot("+", ev.x, ev.y, size=0.5, ctype=ds9.CYAN)

                continue
            
            grd = table.grd
            print "RHL", grd
            assert grd in grades

            if display:
                ds9.line([(ev.x - size, ev.y - size),
                          (ev.x + size, ev.y - size),
                          (ev.x + size, ev.y + size),
                          (ev.x - size, ev.y + size),
                          (ev.x - size, ev.y - size)], frame=0, ctype=ctypes[grd])
                    
    #table.dump_head()
    if False:
        table.dump_hist()
        
    if False:
        table.dump_table()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

thresh=30
fileName="/Users/rhl/TeX/Talks/LSST/Camera-2012/HandsOn/Fe55/Data/C0_20090717-214738-141.fits.gz"

if __name__ == "__main__":
    processImage(thresh=thresh, fileName=fileName)
