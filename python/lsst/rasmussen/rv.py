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
import numpy as np
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
    showMask = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def processImage(thresh, fileName, grades=range(9), searchThresh=None, split=None, emulateMedpict=False,
                 outputFile=None,
                 showRejects=False, showUnknown=False, showGrades=True):
    image = afwImage.ImageF(fileName)

    image -= afwMath.makeStatistics(image, afwMath.MEDIAN).getValue()

    if emulateMedpict:
        sim = image.Factory(image, afwGeom.BoxI(afwGeom.PointI(2008, 0), image.getBBox().getMax()))
        sim.set(0)

    if searchThresh is None:
        searchThresh = thresh

    fs = afwDetect.FootprintSet(image, afwDetect.Threshold(searchThresh))

    if display:
        if showMask:
            mi = afwImage.makeMaskedImage(image)
            afwDetect.setMaskFromFootprintList(mi.getMask(), fs.getFootprints(), 0x4)
        else:
            mi = image
        ds9.mtv(mi, title="bkgd subtracted", frame=0)

    events = []
    for foot in fs.getFootprints():
        for i, peak in enumerate(foot.getPeaks()):
            x, y = peak.getIx(), peak.getIy()

            if x >= 2008:               # ignore overclock
                continue

            if False:
                if x == 793 and y == 59:
                    ds9.pan(x, y)
                    import pdb; pdb.set_trace() 
            
            if i > 0:
                if np.hypot(x - x0, y - y0) < 0: # XXX
                    x0, y0 = x, y
                    continue

            x0, y0 = x, y
            #
            # medpict fails to find some events with two identical adjacent peak values,
            # in particular it uses this peak criterion
            #
            # If you make emulateMedpict == False and run showMedpict(dmEvents)
            # you'll see the real events it missed
            #
            if emulateMedpict:
                v00 = image.get(x, y)
                try:
                    if not (v00 >  image.get(x - 1, y - 1) and
                            v00 >  image.get(x    , y - 1) and
                            v00 >  image.get(x + 1, y - 1) and
                            v00 >  image.get(x - 1, y    ) and
                            v00 >= image.get(x + 1, y    ) and
                            v00 >= image.get(x - 1, y + 1) and
                            v00 >= image.get(x    , y + 1) and
                            v00 >= image.get(x + 1, y + 1)):
                        continue
                except lsst.pex.exceptions.LsstCppException, e:
                    pass

            try:
                events.append(ras.Event(image, afwGeom.PointI(x, y)))
            except lsst.pex.exceptions.LsstCppException, e:
                pass

    print "Found %d events" % (len(events))

    if split is None:
        split = int(0.33*thresh)
    filt = sum([1 << g for g in grades])
    #
    # Pass through Gflt filter
    #
    table = ras.HistogramTableGflt(filt, thresh, split)

    status = []
    for i, ev in enumerate(events):
        success = table.process_event(ev)
        grd = table.grd
        status.append([success, grd])

        if success:
            print "%d %d %d %d %g %d" % (ev.x, ev.y, table.grd, table.sum, ev[4], 0 if True else table.p9)
            
    #
    # and Xygpx filter
    #
    table = ras.HistogramTableXygpx(ras.HistogramTableXygpx.p_9, filt, thresh, split)

    for i, ev in enumerate(events):
        success = table.process_event(ev)
        grd = table.grd
        status[i].append(success)

    size = 1.6                          # half-size of box to draw
    ctypes = dict(zip(range(-1, 8), (ds9.WHITE,
                                     ds9.RED, ds9.GREEN, ds9.BLUE,      # Grade: 0 1 2
                                     ds9.CYAN, ds9.MAGENTA, ds9.YELLOW, #        3 4 5
                                     ds9.WHITE, ds9.RED)))              #        6 7

    with ds9.Buffering():
        for ev, st in zip(events, status):
            successGflt, grd, successXygpx = st
            success = successXygpx

            if not success:
                if display and (grd == -1 and showUnknown or grd >= 0 and showRejects):
                    ds9.dot("+", ev.x, ev.y, size=0.5, ctype=ctypes[grd])
                    if showGrades:
                        ds9.dot(str(grd), ev.x + 1.5, ev.y, frame=0, ctype=ctypes[grd])
                    
                continue

            #assert grd in grades

            if display:
                ds9.line([(ev.x - size, ev.y - size),
                          (ev.x + size, ev.y - size),
                          (ev.x + size, ev.y + size),
                          (ev.x - size, ev.y + size),
                          (ev.x - size, ev.y - size)], frame=0, ctype=ctypes[grd])
                if showGrades:
                    ds9.dot(str(grd), ev.x + size + 1, ev.y - size, frame=0, ctype=ctypes[grd])

    with open("XXX", "w") as fd:
        table.dump_head(fd)
        table.dump_hist(fd)
    #table.dump_table()

    return table, image, events

def showMedpict(fileName="events.dat", events=None, image=None):
    medpictEvents = ras.readEventFile(fileName)

    #
    # Look for events that medpict missed
    #
    if events:
        x = np.empty(len(medpictEvents))
        y = np.empty(len(medpictEvents))

        for i, ev in enumerate(medpictEvents):
            x[i] = ev.x
            y[i] = ev.y

        for ev in events:
            if events:
                d = np.hypot(x - ev.x, y - ev.y)
                dmin = np.min(d)
                if dmin > 0:
                    print "medpict missed:", " ".join([str(_) for _ in zip(x[d == dmin], y[d == dmin])])

                    ds9.dot("o", ev.x, ev.y, size=5, ctype=ds9.GREEN)
                    if False:
                        ds9.pan(ev.x, ev.y)
                        ds9.flush()
                        import pdb; pdb.set_trace() 
        #
        # Look for events that the DM stack missed
        #
        if events:
            x = np.empty(len(events))
            y = np.empty(len(events))

            for i, ev in enumerate(events):
                x[i] = ev.x
                y[i] = ev.y

    with ds9.Buffering():
        for ev in medpictEvents:
            if events:
                d = np.hypot(x - ev.x, y - ev.y)
                dmin = np.min(d)
                if dmin > 0:
                    print "DM missed:", " ".join([str(_) for _ in zip(x[d == dmin], y[d == dmin])])
                    
                    ds9.dot("o", ev.x, ev.y, size=5, ctype=ds9.RED)
                    if False:
                        ds9.pan(ev.x, ev.y)
                        ds9.flush()
                        import pdb; pdb.set_trace() 

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

searchThresh=20
thresh=30
fileName="/Users/rhl/TeX/Talks/LSST/Camera-2012/HandsOn/Fe55/Data/C0_20090717-214738-141.fits.gz"

if __name__ == "__main__":
    processImage(searchThresh=searchThresh, thresh=thresh, fileName=fileName)

    
