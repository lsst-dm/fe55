import numpy as np
import lsst.daf.base as dafBase
import lsst.pex.exceptions as pexExcept
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9

def makeAmp(md, channelNo=None, trim=True):
    """Make a cameraGeom.Amp from metadata (or ab initio given it's channel number)"""
    #
    # Numbers that aren't in the header
    #
    nRow = 2000                         # number of real rows of CCD pixels readout through an amp
    nCol = 512                          # number of real columns of CCD pixels
    data0 = 10                          # first column of real pixels
    gain, readNoise, saturationLevel = 1.0, 0, 0

    if md is not None:
        channelNo = md.get("CHANNEL")
        ewidth, eheight = [int(_[2:]) for _ in md.get("DATASEC")[1:-1].split(",")]

        ltv = np.matrix((md.get("LTV1"), md.get("LTV2"))).getT()

        ltm = np.matrix(np.zeros((2, 2)))
        for i in range(2):
            for j in range(2):
                try:
                    ltm[i, j] = md.get("LTM%d_%d" % (j + 1, i + 1), 0.0)
                except pexExcept.LsstException:
                    pass

        if abs(ltm[0, 0]) != 1.0 or abs(ltm[1, 1]) != 1.0:
            raise RuntimeError("I refuse to handle non-square pixels: LTM = %s" % ltm)        
        if ltm[1, 0] != 0.0 or ltm[0, 1] != 0.0:
            raise RuntimeError("I refuse to handle sheared detectors: LTM = %s" % ltm)

        rotate90 = 0                    # no need to rotate the amp to fit in the CCD
        flipLR = False                  # no need to flip the amp horizontally
        if ltm[0, 0] < 0:
            if ltm[1, 1] < 0:
                rotate90 = 2
                flipLR = True
        else:
            if ltm[1, 1] < 0:
                rotate90 = 2
                flipLR = True
        #
        # I think the camera team got the geometry wrong
        # For channel 1 they have:
        #   DATASEC : [1:542,1:2022]
        #   DETSEC  : [0:542,0:2022]  (n.b. that's  543x2023)
        #   DETSIZE : [0:4336,0:4044] (n.b. that's 4337x4045)
        # If the 0 is correct, it implies that there's logically 1 pre-scan pixel so Cx starts at 0.0
        # Note also that channel 16 has
        #   LTV1 :  542
        #   LTV2 : 4044
        # which means that the pixel with (Ic, Il) == (542, 4044) maps to (Cx, Cy) = (0.0, 0.0)
        # The problems with this theory are:
        #   1. There are more than 1 overscan pixels (I think 10)
        #   2. The size of the DETSEC and DETSIZE are wrong --- they should be e.g. [0:541, 0:2021]
        #
        # I think a more likely problem is that there's an off-by-one error in LTV[12] for the rotated
        # channels at the top of the chip; I shall fix things with this assumption
        #
        if rotate90 == 2:
            ltv += 1

        iltm = ltm.getI()
        llc = np.matrix((1.0, 1.0)).getT()
        urc = np.matrix((ewidth, eheight)).getT()
        llc, urc = iltm*(llc - ltv), iltm*(urc - ltv)

        if rotate90 == 2:
            llc, urc = urc, llc

        llc -= 1; urc -= 1                  # convert to C/C++/python/LSST 0-indexed convention

        iCol, iRow = int(llc[0]//ewidth), int(llc[1]//eheight) # int() converts from numpy.matrix
    else:
        ewidth, eheight = 542, 2022

        if channelNo <= 8:
            iCol, iRow = (channelNo - 1), 0
            rotate90 = 0                # Amp image is in the
            flipLR = False              #    same orientation as the CCD image
        elif channelNo <= 16:
            iCol, iRow = (16 - channelNo), 1

            rotate90 = 2                # Amp image is rotated
            flipLR = True               #     and flipped left-right relative to the CCD image
        else:
            raise StopIteration()

    dataSec = afwGeom.BoxI(afwGeom.PointI(data0, 0), afwGeom.ExtentI(nCol, nRow))
    biasSec = afwGeom.BoxI(afwGeom.PointI(data0 + nCol, 0), afwGeom.PointI(ewidth - 1, nRow - 1))
    biasSec = afwGeom.BoxI(afwGeom.PointI(data0 + nCol, 0), afwGeom.PointI(ewidth - 1, nRow - 1))
    allPixelsInAmp = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.ExtentI(ewidth, eheight))

    eParams = cameraGeom.ElectronicParams(gain, readNoise, saturationLevel)

    amp = cameraGeom.Amp(cameraGeom.Id(channelNo, "ID%d" % channelNo, iCol, iRow),
                         allPixelsInAmp, biasSec, dataSec, eParams)
    #The following maps how the amp pixels must be changed to go from electronic (on disk) coordinates
    #to detector coordinates.  This also sets the readout corner.
    amp.setElectronicToChipLayout(afwGeom.Point2I(iCol, iRow), rotate90, flipLR, cameraGeom.Amp.AMP)

    amp.setTrimmed(trim)
    return amp

def makeCcd(fileName, fromHeader=False, trim=True):
    ccd = cameraGeom.Ccd(cameraGeom.Id(0))

    a = 0                               # one-less than the next HDU
    while True:                         # while there are valid HDUs
        a += 1
        if fromHeader:
            try:
                hdu = 1 + a
                md, channelNo = afwImage.readMetadata(fileName, hdu), None
            except lsst.pex.exceptions.LsstCppException:
                if hdu == 1:            # an empty PDU
                    continue
                break
            
            if "TTYPE1" in md.names():  # not an image
                break
        else:
            md, channelNo = None, a
        #
        # Add amp to the Ccd
        #
        try:
            ccd.addAmp(makeAmp(md, channelNo, trim))
        except StopIteration:
            break

    return ccd

def assembleCcd(fileName, trim=False, perRow=True):
    """Assemble a complete CCD image.  If trim is true, bias subtract and trim to the "real" pixels
If perRow is True, estimate the bias level for each row of the overclock

Return a tuplle of (afwCameraGeom.Ccd, ccdImage)
    """
    ccd = makeCcd(fileName, trim=trim)

    ccdImage = afwImage.ImageF(ccd.getAllPixels(trim))

    for a in ccd:
        channelNo = a.getId().getSerial()
        im = afwImage.ImageF(fileName, 1 + channelNo)

        if trim:
            bias = im.Factory(im, a.getDiskBiasSec())
            im = im.Factory(im, a.getDiskDataSec())

            if perRow:
                data = im.getArray();
                biasVec = np.median(bias.getArray(), True)
                for i in range(data.shape[1]):
                    data[:, i] -= biasVec
            else:
                im -= afwMath.makeStatistics(bias, afwMath.MEANCLIP).getValue()

        a.setTrimmed(True)

        sub = ccdImage.Factory(ccdImage, a.getAllPixels(trim))
        sub <<= a.prepareAmpData(im)

    return ccd, ccdImage
    
if __name__ == "__main__":
    makeCcd("/Users/rhl/Desktop/fe55_0600s.fits.gz")
