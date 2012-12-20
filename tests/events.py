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
Tests for Event

Run with:
   ./events.py
or
   python
   >>> import events; events.run()
"""

import os
import os.path

import sys
import unittest

import lsst.utils.tests as utilsTests
import lsst.pex.exceptions
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9

import lsst.rasmussen as ras

import numpy

try:
    type(display)
except NameError:
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class EventsTestCase(unittest.TestCase):
    """A test case for Image"""
    def setUp(self):
        self.image = afwImage.ImageF(afwGeom.ExtentI(100, 200))
        self.centers = [afwGeom.PointI(10, 30),
                        afwGeom.PointI(50, 100),
                        ]

        self.table = ras.HistogramTableGflt()

    def tearDown(self):
        del self.image
        del self.centers
        del self.table

    def testCtor(self):
        events = [ras.Event(self.image, ev) for ev in self.centers]

        ev = events[0]

        def offChip():
            ras.Event(self.image, self.image.getBBox().getMax() + afwGeom.ExtentI(10, 10))

        utilsTests.assertRaisesLsstCpp(self, lsst.pex.exceptions.OutOfRangeException, offChip)

    if False:
        def testEventTable_dump_table(self):
            self.table.dump_table()

        def testEventTable_dump_head(self):
            self.table.dump_head()

        def testEventTable_dump_hist(self):
            self.table.dump_hist()
    else:
        print >> sys.stderr, "Skipping table.dump* tests"
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(EventsTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)