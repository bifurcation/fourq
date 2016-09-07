#!/usr/bin/env python

import re
from fields import GFp2

def fmtpt(P):
    fmt = "(0x{:032x} + 0x{:032x} * i, 0x{:032x} + 9x{:032x} * i)"
    return fmt.format(P[0][0], P[0][1], P[1][0], P[1][1])

# Take a curve in R1 or R4 and clear the projective factor
def toAffine(P):
    if len(P) == 2:
        return P
    if len(P) != 3 and len(P) != 5:
        raise Exception("Representation unsupported for normalization")
    zi = GFp2.inv(P[2])
    xz = GFp2.mul(P[0], zi)
    yz = GFp2.mul(P[1], zi)
    return (xz, yz)

def test(label, sample, ref):
    if sample == ref:
        print "[PASS] {}".format(label)
    else:
        print "[FAIL] {} {}".format(label, sample)

def testpt(label, sample, ref):
    aff1 = toAffine(sample)
    aff2 = toAffine(ref)
    if aff1 == aff2:
        print "[PASS] {}".format(label)
    else:
        print "[FAIL] {} {}".format(label, fmtpt(aff1))
