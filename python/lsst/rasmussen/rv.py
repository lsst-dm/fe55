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
import lsst.daf.base as dafBase
import lsst.pex.exceptions
import lsst.afw.cameraGeom as afwCameraGeom
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9

import lsst.rasmussen as ras

import numpy

try:
    import matplotlib.pyplot as plt
    try:
        fig
    except NameError:
        fig = None
except ImportError:
    plt = None

try:
    type(display)
except NameError:
    ctypes = dict(zip(range(-1, 8),
                      (ds9.WHITE,   # UNKNOWN
                       ds9.RED,     # 0 SINGLE                 Single pixel
                       ds9.GREEN,   # 1 SINGLE_P_CORNER        Single pixel + corner(s)
                       ds9.BLUE,    # 2 VERTICAL_SPLIT         Vertical split (+ detached corner(s))
                       ds9.CYAN,    # 3 LEFT_SPLIT             Left split (+ detached corner(s))
                       ds9.MAGENTA, # 4 RIGHT_SPLIT            Right split (+ detached corner(s))
                       ds9.YELLOW,  # 5 SINGLE_SIDED_P_CORNER  Single-sided split + touched corner
                       ds9.WHITE,   # 6 ELL_SQUARE_P_CORNER    L or square (+ detached corner)
                       ds9.RED      # 7 OTHER                  all others
                       )))

    display = False
    showMask = False

def makeAmp(md, emulateMedpict=False):
    startData      = 0    if emulateMedpict else 50
    startOverclock = 2008 if emulateMedpict else 2065

    gain, readNoise, saturation = 1.0, 5.0, 40000
    eparams = afwCameraGeom.ElectronicParams(gain, readNoise, saturation)

    allPixels = afwGeom.BoxI(afwGeom.Point2I(0, 0), afwGeom.ExtentI(2116, 1000))
    biasSec = afwGeom.BoxI(afwGeom.PointI(startOverclock, 0),
                           allPixels.getMax())
    dataSec = afwGeom.BoxI(afwGeom.PointI(startData, 0),
                           afwGeom.PointI(startOverclock - 1, allPixels.getMaxY()))
    #
    # Check those values from the header... except that they are absent so we can't
    #
    assert allPixels.getWidth() == md.get("WIDTH")
    assert allPixels.getHeight() == md.get("HEIGHT")
    ccdSer = md.get("CCD_SER") if False else 0 # not set

    return afwCameraGeom.Amp(afwCameraGeom.Id(ccdSer), allPixels, biasSec, dataSec, eparams)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def processImage(thresh, fileName, grades=range(8), searchThresh=None, split=None,
                 emulateMedpict=False, outputHistFile=None, outputEventsFile=None,
                 showRejects=False, showUnknown=False, showGrades=True, plot=True, subplots=False):

    md = dafBase.PropertyList()
    image = afwImage.ImageF(fileName, 1, md)
    amp = makeAmp(md)

    bias = image.Factory(image, amp.getBiasSec())

    if emulateMedpict:
        image -= afwMath.makeStatistics(image, afwMath.MEDIAN).getValue()
        bias.set(0)
    else:
        image -= afwMath.makeStatistics(bias, afwMath.MEDIAN).getValue()

    if searchThresh is None:
        searchThresh = thresh

    dataSec = image.Factory(image, amp.getDataSec())
    fs = afwDetect.FootprintSet(dataSec, afwDetect.Threshold(searchThresh))

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
            #
            # medpict fails to find some events with two identical adjacent peak values,
            # in particular it uses th peak criterion embodied in the following logic
            #
            # If you set emulateMedpict == False and run showMedpict(dmEvents)
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

    print "Found  %5d events" % (len(events))

    if split is None:
        split = int(0.33*thresh)
    filt = sum([1 << g for g in grades])
    reset = 0.0
    table = ras.HistogramTableBase(thresh, split,
                                   ras.HistogramTableBase.T1, reset)
    table.setFilter(filt)
    table.setCalctype(calcType)

    status = []
    for i, ev in enumerate(events):
        success = table.process_event(ev)
        status.append([success, table.sum])
    print "Passed %5d events" % (sum([_[0] for _ in status]))

    if outputEventsFile:
        with open(outputEventsFile, "w") as fd:
            for stat, ev in zip(status, events):
                stat, _sum = stat
                if stat:
                    print >> fd, "%d %d %d %d %g %d" % (ev.x, ev.y, ev.grade, _sum, ev[4], ev.p9)

    size = 1.6                          # half-size of box to draw

    with ds9.Buffering():
        for ev, st in zip(events, status):
            success, _sum = st

            if not success:
                if display and (ev.grade == -1 and showUnknown or ev.grade >= 0 and showRejects):
                    ds9.dot("+", ev.x, ev.y, size=0.5, ctype=ctypes[ev.grade])
                    if showGrades:
                        ds9.dot(str(ev.grade), ev.x + 1.5, ev.y, frame=0, ctype=ctypes[ev.grade])
                    
                continue

            if display:
                ds9.line([(ev.x - size, ev.y - size),
                          (ev.x + size, ev.y - size),
                          (ev.x + size, ev.y + size),
                          (ev.x - size, ev.y + size),
                          (ev.x - size, ev.y - size)], frame=0, ctype=ctypes[ev.grade])
                if showGrades:
                    ds9.dot(str(ev.grade), ev.x + size + 1, ev.y - size, frame=0, ctype=ctypes[ev.grade])

    if outputHistFile:
        with open(outputHistFile, "w") as fd:
            table.dump_head(fd)
            table.dump_hist(fd)
    #table.dump_table()

    if plot:
        plot_hist(table, title="Event = %g Split = %g Source = %s" % (thresh, split, "unknown"),
                  subplots=subplots)

    return table, image, events

def plot_hist(table, title=None, subplots=False, scaleSubplots=False):
    """Plot the histograms in a HistogramTable

    \param table The table to plot
    \param title A title for the plot, if not None
    \param subplots Make a separate plot for each event grade
    \param scaleSubplots Scale each subplot separately
    """
    if not plt:
        return
    
    global fig
    if fig is None:
        fig = plt.figure()
    else:
        fig.clf()
    if subplots:
        plt.subplots_adjust(0.1, 0.1, 0.95, 0.95, wspace=0.0, hspace=0.0)
    else:
        axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    if subplots:
        nplot = sum([sum(_) != 0 for _ in table.histo])
    else:
        nplot = 1

    x = np.arange(0, table.MAXADU)
    i = 0
    yMax = 0
    for g, label in enumerate(["N(S)",
                               "N(S+)",
                               "N(Pv)",
                               "N(Pl)",
                               "N(Pr)",
                               "N(P+)",
                               "N(L+Q)",
                               "N(O)",]):
        y = table.histo[g]
        if sum(y) == 0:
            continue

        dy = np.sqrt(y)                 # error in y

        color="black" if ctypes[g] == "white" else ctypes[g]

        i += 1
        if subplots:
            axes = fig.add_subplot(nplot, 1, i)
            
        axes.set_yscale('log', nonposy='clip')

        lineLabel = "%d %s" % (g, label)
        axes.step(x, y, label=lineLabel, where='mid', color=color)
        #
        # Draw an error band.  The "incomprehensible list comprehensions" flatten the zips:
        #    [_ for foo in goo for _ in foo]
        # reads as:
        #    for foo in goo:
        #       for _ in foo:
        #          _
        # allowing the band to follow the histogram, not its centres; xx = -0.5, 0.5, 0.5, 1.5, 1.5 ...
        xx =     np.array([_ for xmp in zip(x - 0.5, x + 0.5) for _ in xmp])
        ym =     np.array([_ for _y in y - dy for _ in (_y, _y)])
        yp =     np.array([_ for _y in y + dy for _ in (_y, _y)])

        alpha = 0.35
        axes.fill_between(xx, ym, yp, color=color, alpha=alpha)

        if subplots:
            if i < nplot:
                axes.set_xticklabels([])
            axes.text(0.85, 0.70, lineLabel, ha="left", transform=axes.transAxes)
        else:
            axes.set_xlabel("Pulse Height (ADU)")

        axes.set_ylabel("lg(N)")
        axes.set_xlim(table.min_adu - 10, table.max_adu + 10)
        ytop = axes.get_ylim()[1]
        if ytop > yMax:
            yMax = ytop

    if subplots:
        if scaleSubplots:
            yMax = None
        else:
            yMax *= 0.99                    # Labels at the top of the y-axis overlap
        for i in range(nplot):
            fig.add_subplot(nplot, 1, i + 1).set_ylim(0.8, yMax)
    else:
        axes.set_ylim(0.8, yMax)

    if title:
        fig.suptitle(title)
        
    if not subplots:
        axes.legend(loc=1)
    fig.show()

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

    
