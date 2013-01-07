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
import lsst.rasmussen.cameraGeom as cameraGeom

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
    type(ctypes)
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

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def processImage(thresh, fileNames, grades=range(8), searchThresh=None, split=None,
                 calcType=ras.HistogramTable.P_9,
                 outputHistFile=None, outputEventsFile=None, assembleCcd=False, plotByAmp=False,
                 plot=True, subplots=False, xlim=[None, 650], ylim=[None, None],
                 displayRejects=False, displayUnknown=False, displayGrades=True, display=False, 
                 emulateMedpict=None    # not used
                 ):

    if searchThresh is None:
        searchThresh = thresh

    nImage = 0                          # number of images we've processed
    events = []
    ampIds = set()
    for frameNum, fileName in enumerate(fileNames):
        # Read file
        hdu = 0                         # one-less than the next HDU
        while True:                     # while there are valid HDUs
            hdu += 1

            if assembleCcd:
                if hdu > 1:
                    break
                ccd, image = cameraGeom.assembleCcd(fileName, trim=True)
                dataSec = image
                ampIds = set(_.getId().getSerial() for _ in ccd)
            else:
                ccd = None              # we don't have an assembled Ccd
                md = dafBase.PropertyList()
                try:
                    image = afwImage.ImageF(fileName, hdu, md)
                except lsst.pex.exceptions.LsstCppException:
                    if hdu == 1:            # an empty PDU
                        continue
                    break

                # Get the image's camera geometry (e.g. where is the datasec?)
                amp = cameraGeom.makeAmp(md)
                ampIds.add(amp.getId().getSerial())
                
                # Subtract the bias, estimated as the median of the biassec
                bias = image.Factory(image, amp.getDiskBiasSec())
                image -= afwMath.makeStatistics(bias, afwMath.MEDIAN).getValue()
                # Search the datasec for Fe55 events
                dataSec = image.Factory(image, amp.getDiskDataSec())

            nImage += 1
            fs = afwDetect.FootprintSet(dataSec, afwDetect.Threshold(searchThresh))

            if display:
                if True:
                    ds9.erase()
                else:
                    mi = afwImage.makeMaskedImage(image)
                    afwDetect.setMaskFromFootprintList(mi.getMask(), fs.getFootprints(), 0x4)
                    ds9.mtv(mi, title="bkgd subtracted", frame=0)
                    del mi

            # Convert all the peaks within the detections to Events
            for foot in fs.getFootprints():
                for i, peak in enumerate(foot.getPeaks()):
                    peakPos = afwGeom.PointI(peak.getIx(), peak.getIy())
                    if ccd:
                        amp = ccd.findAmp(peakPos, True)
                        
                    try:
                        events.append(ras.Event(image, peakPos, frameNum, amp.getId().getSerial()))
                    except lsst.pex.exceptions.LsstCppException, e:
                        pass
    #
    # Prepare to go through all our events, building our histograms
    #
    if split is None:
        split = int(0.33*thresh)

    filt = sum([1 << g for g in grades])

    tables = {}
    if plotByAmp:
        for aid in ampIds:
            tables[aid] = ras.HistogramTable(thresh, split)
    else:
        table = ras.HistogramTable(thresh, split)
        for aid in ampIds:
            tables[aid] = table

    for table in tables.values():
        table.setFilter(filt)
        table.setCalctype(calcType)
        table.setReset(ras.HistogramTable.T1, 0.0)
    del table

    # Process the events
    table0 = tables.values()[0]
    if plotByAmp:
        status = len(events)*[None]
        for i, ev in enumerate(events):
            status[i] = tables[ev.chipnum].process_event(ev)
    else:
        status = [table0.process_event(ev) for ev in events]
    print "Passed %5d events" % (sum(status))
    #
    # Done.  Output...
    #
    if outputEventsFile:
        with open(outputEventsFile, "w") as fd:
            for stat, ev in zip(status, events):
                if stat:
                    print >> fd, "%d %d %d %d %g %d" % (ev.x, ev.y, ev.grade, ev.sum, ev[4], ev.p9)

    if outputHistFile:
        with open(outputHistFile, "w") as fd:
            table0.dump_head(fd, "unknown", sum(status))
            table0.dump_hist(fd)
    #table0.dump_table()

    if nImage > 1:
        if display:
            print >> sys.stderr, "No point plotting events as they refer to multiple amps/ccds"
    elif display:
        size = 1.6                      # half-size of box to draw
        with ds9.Buffering():
            for ev, success in zip(events, status):
                if not success:
                    if ev.grade == -1 and displayUnknown or ev.grade >= 0 and displayRejects:
                        ds9.dot("+", ev.x, ev.y, size=0.5, ctype=ctypes[ev.grade])
                        if displayGrades:
                            ds9.dot(str(ev.grade), ev.x + 1.5, ev.y, frame=0, ctype=ctypes[ev.grade])
                    
                    continue

                ds9.line([(ev.x - size, ev.y - size),
                          (ev.x + size, ev.y - size),
                          (ev.x + size, ev.y + size),
                          (ev.x - size, ev.y + size),
                          (ev.x - size, ev.y - size)], frame=0, ctype=ctypes[ev.grade])
                if displayGrades:
                    ds9.dot(str(ev.grade), ev.x + size + 1, ev.y - size, frame=0, ctype=ctypes[ev.grade])

    if plot:
        plot_hist(tables, title="Event = %g Split = %g Source = %s N=%d" %
                  (thresh, split, os.path.basename(fileName), sum(status)),
                  xlim=xlim, ylim=ylim, subplots=subplots)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plot_hist(tables, title=None, subplots=False, scaleSubplots=False, xlim=[None, None], ylim=[None, None]):
    """Plot the histograms in a HistogramTable

    \param table The table to plot
    \param title A title for the plot, if not None
    \param subplots Make a separate plot for each event grade
    \param scaleSubplots Scale each subplot separately
    """
    if not plt:
        return

    nTable = len(set(tables.values()))
    table0 = tables.values()[0]
    if nTable == 1:
        tables = {None: table0}
        
    if nTable > 1:         # more than just a single table
        if subplots:
            print >> sys.stderr, "You may not request per-grade subplots and also provide per-amp histograms"
        subplots = True                 # but now it means per-amp, not per grade

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
        if nTable == 1:
            nPlot = sum([sum(_) != 0 for _ in table0.histo])
        else:
            nPlot = nTable
    else:
        nPlot = 1

    if False:
        nRow = 6
        nCol = nPlot//nRow
        while nRow*nCol < nPlot:
            nCol += 1
    else:
        nCol = 2
        nRow = nPlot//nCol
        while nRow*nCol < nPlot:
            nRow += 1

    x = np.arange(0, table0.MAXADU)
    i = 0
    xMax, yMax = 0, 0
    for tableKey in sorted(tables.keys()):
        table = tables[tableKey]
        newPanel = True
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

            if newPanel:
                i += 1
                if subplots:
                    axes = fig.add_subplot(nRow, nCol, i)

            axes.set_yscale('log', nonposy='clip')

            if nTable == 1:
                lineLabel = "%d %s" % (g, label)
            else:
                lineLabel = "Amp %d" % (tableKey)

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
                if i + nCol <= nPlot:
                    axes.set_xticklabels([])
                if newPanel:
                    axes.text(0.70, 0.70, lineLabel, ha="left", transform=axes.transAxes)
            else:
                axes.set_xlabel("Pulse Height (ADU)")

            if i%nCol == 1:
                axes.set_ylabel("lg(N)")
            else:
                axes.set_yticklabels([])
                
            axes.set_xlim(table.min_adu - 10, table.max_adu + 10)
            xtop, ytop = axes.get_xlim()[1], axes.get_ylim()[1]
            if xtop > xMax:
                xMax = xtop
            if ytop > yMax:
                yMax = ytop

            if nTable > 0:
                newPanel = False        # don't get a new panel for every grade

    if scaleSubplots:
        yMax = None
    elif subplots:
        yMax *= 0.99                    # Labels at the top of the y-axis overlap


    if xlim[0] is None: xlim[0] = 0
    if xlim[1] is None: xlim[1] = xMax
    if ylim[0] is None: ylim[0] = 0.8
    if ylim[1] is None: ylim[1] = yMax
    if subplots:
        for i in range(nPlot):
            fig.add_subplot(nRow, nCol, i + 1).set_xlim(*xlim)
            fig.add_subplot(nRow, nCol, i + 1).set_ylim(*ylim)
    else:
        axes.set_xlim(*xlim)
        axes.set_ylim(*ylim)

    if title:
        fig.suptitle(title)
        
    if not subplots:
        axes.legend(loc=1)
    fig.show()

def calcTypeFromString(string):
    P_names =  [_ for _ in dir(ras.HistogramTable) if _.startswith("P_")]
    if string not in P_names:
        raise RuntimeError("Name %s is not in [%s]" % (string, ", ".join(P_names)))

    return eval("ras.HistogramTable.%s" % string)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

searchThresh=20
thresh=30
fileName="/Users/rhl/TeX/Talks/LSST/Camera-2012/HandsOn/Fe55/Data/C0_20090717-214738-141.fits.gz"

if __name__ == "__main__":
    processImage(searchThresh=searchThresh, thresh=thresh, fileName=[fileName])
