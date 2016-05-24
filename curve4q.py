#!/usr/bin/env python

from random import getrandbits
from time import time
import re

# For convenience
zero = (0, 0)
one  = (1, 0)

# For conditional swaps (128-bit masks)
CFALSE = 0x00000000000000000000000000000000
CTRUE  = 0xffffffffffffffffffffffffffffffff

# Field characteristic
p = 0x7fffffffffffffffffffffffffffffff

# Curve parameter as a field element tuple
dx = 0x5e472f846657e0fcb3821488f1fc0c8dL
dy = 0xe40000000000000142L
d = (dy, dx)

# Generator from FourQ.c, in extended twisted Edwards coordinates
Gx = (0x1A3472237C2FB305286592AD7B3833AA, 0x1E1F553F2878AA9C96869FB360AC77F6)
Gy = (0x0E3FEE9BA120785AB924A2462BCBB287, 0x6E1C4AF8630E024249A7C344844C8B5C)
G = (Gx, Gy, one, Gx, Gy)

# Neutral / zero element of the curve
O = (zero, one, one, zero, one)


########## Field Arithmetic ##########

# GF(p^2) elements are represented as 2-tuples of longs
#  x + y i == (x, y)
#
# GF(p) operations are just done using python longs, except for inversion:
#
# a^-1 = a^{p-2}
#
# (a0 + a1 i) + (b0 + b1 i) = (a0 + b0) + (a1 + b1) i
# (a0 + a1 i) * (b0 + b1 i) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0) i
# (a0 + a1 i)^2 = (a0^2 - a1^2) + (2 a0 a1) i
# (a0 + a1 i)^-1 = (a0 - a1 i) / (a0^2 + a1^2)

def fpadd(a, b):
    return (a + b) % p

def fpsub(a, b):
    return (a - b) % p

def fpmul(a, b):
    return (a * b) % p

def fpsqr(a):
    return (a * a) % p

def fpinv(a):
    t2 = fpsqr(a)           # 2
    t2 = fpmul(a, t2)       # 3
    t3 = fpsqr(t2)          # 6
    t3 = fpsqr(t3)          # 12
    t3 = fpmul(t2, t3)      # 15
    t4 = fpsqr(t3)          # 30
    t4 = fpsqr(t4)          # 60
    t4 = fpsqr(t4)          # 120
    t4 = fpsqr(t4)          # 240
    t4 = fpmul(t3, t4)      # 2^8 - 2^0
    t5 = fpsqr(t4)          # 2^9 - 2^1
    for i in range(7):      # 2^16 - 2^8
        t5 = fpsqr(t5)
    t5 = fpmul(t4, t5)      # 2^16 - 2^0
    t2 = fpsqr(t5)          # 2^17 - 2^1
    for i in range(15):     # 2^32 - 2^16
        t2 = fpsqr(t2)
    t2 = fpmul(t5, t2)      # 2^32 - 2^0
    t1 = fpsqr(t2)          # 2^33 - 2^1
    for i in range(31):     # 2^64 - 2^32
        t1 = fpsqr(t1)
    t1 = fpmul(t1, t2)      # 2^64 - 2^0
    for i in  range(32):    # 2^96 - 2^32
        t1 = fpsqr(t1)
    t1 = fpmul(t1, t2)      # 2^96 - 2^0
    for i in  range(16):    # 2^112 - 2^16
        t1 = fpsqr(t1)
    t1 = fpmul(t1, t5)      # 2^112 - 2^0
    for i in  range(8):     # 2^120 - 2^8
        t1 = fpsqr(t1)
    t1 = fpmul(t1, t4)      # 2^120 - 2^0
    for i in  range(4):     # 2^124 - 2^4
        t1 = fpsqr(t1)
    t1 = fpmul(t1, t3)      # 2^124 - 2^0
    t1 = fpsqr(t1)          # 2^125 - 2^1
    t1 = fpmul(t1, a)       # 2^125 - 2^0
    t1 = fpsqr(t1)          # 2^126 - 2^1
    t1 = fpsqr(t1)          # 2^127 - 2^2
    return fpmul(t1, a)     # 2^127 - 3

def fpneg(a):
    return (p - a) % p

def fpcswap(c, x, y):
    return (((x ^ y) & c) ^ y)

###

# Counters for GF(p^2) operations

class GFp2:
    A = 0
    S = 0
    M = 0
    I = 0

    @staticmethod
    def ctr_reset():
        GFp2.A = 0
        GFp2.S = 0
        GFp2.M = 0
        GFp2.I = 0

    @staticmethod
    def ctr():
        return (GFp2.A, GFp2.S, GFp2.M, GFp2.I)

    @staticmethod
    def add(a, b):
        GFp2.A += 1
        return (fpadd(a[0], b[0]), fpadd(a[1], b[1]))

    @staticmethod
    def sub(a, b):
        GFp2.A += 1
        return (fpsub(a[0], b[0]), fpsub(a[1], b[1]))

    @staticmethod
    def addsub(a, b):
        return GFp2.sub(GFp2.add(a, a), b)

    @staticmethod
    def mul(a, b):
        GFp2.M += 1
        a0b0 = fpmul(a[0], b[0])
        a1b0 = fpmul(a[1], b[0])
        a0b1 = fpmul(a[0], b[1])
        a1b1 = fpmul(a[1], b[1])
        return (fpsub(a0b0, a1b1), fpadd(a0b1, a1b0))

    @staticmethod
    def sqr(a):
        GFp2.S += 1
        a02   = fpmul(a[0], a[0])
        a12   = fpmul(a[1], a[1])
        a0a12 = fpmul(2, fpmul(a[0], a[1]))
        return (fpsub(a02, a12), a0a12)

    @staticmethod
    def neg(a):
        GFp2.A += 1
        return (fpneg(a[0]), fpneg(a[1]))

    @staticmethod
    def conj(a):
        GFp2.A += 0.5
        return (a[0], fpneg(a[1]))

    @staticmethod
    def inv(a):
        GFp2.I += 1
        invmag = fpinv(fpadd(fpsqr(a[0]), fpsqr(a[1])))
        GFp2.M -= 1 # To avoid double-counting
        GFp2.A -= 1
        return GFp2.mul((invmag, 0), GFp2.conj(a))

    # Constant time
    # switch (c) {
    #   case CTRUE:  return x;
    #   case CFALSE: return y;
    # }
    @staticmethod
    def cswap(c, x, y):
        return (fpcswap(c, x[0], y[0]), fpcswap(c, x[1], y[1]))

    @staticmethod
    def fmt(a):
        return "{:032x}:{:032x}".format(a[0], a[1])

###

def rawhex(x):
    return re.sub(r'L$', "", hex(x)[2:])

def fmtpt(P):
    return " : ".join(map(lambda x: rawhex(x[1]) + rawhex(x[0]), P))

def test(label, sample, ref):
    if sample == ref:
        print "[PASS] {}".format(label)
    else:
        print "[FAIL] {} {}".format(label, sample)

def testpt(label, sample, ref):
    A = normalize(sample)
    B = normalize(ref)
    if A[0] == B[0] and A[1] == B[1]:
        print "[PASS] {}".format(label)
    else:
        print "[FAIL] {} {}".format(label, fmtpt(A))

# XXX Not very systematic, but good enough (?)
def fp2_test():
    i8 = fpinv(13)
    test("inv", fpmul(i8, 13), 1)

    one = (1,0)
    i = (0,1)
    x23 = (2, 3)
    x57 = (5, 7)

    test("1+i", GFp2.add(one,i), (1,1))
    test("1*i", GFp2.mul(one,i), (0,1))

    test("add",     GFp2.add(x23, x57), (7, 10))
    test("sub-pos", GFp2.sub(x57, x23), (3, 4))
    test("sub-neg", GFp2.sub(x23, x57), (p-3, p-4))
    test("mul",     GFp2.mul(x23, x57), (p-11, 29))
    test("sqr",     GFp2.sqr(x23), (p-5, 12))

    x23c = GFp2.conj(x23)
    x23i = GFp2.inv(x23)
    test("conj", x23c, (2, fpneg(3)))
    test("inv", GFp2.mul(x23, x23i), one)

    swap = GFp2.cswap(CTRUE, x23, x57)
    noswap = GFp2.cswap(CFALSE, x23, x57)
    test("swap", swap, x23)
    test("noswap", noswap, x57)

#GFp2._test()

########## Point Representations ##########

# Points are represented as tuples in any of the forms in Table 2
#
# R1    X   Y   Z   Ta  Tb
# R2    X+Y Y-X 2Z  2dT
# R3    X+Y Y-X Z   T
# R4    X   Y   Z
#
# Expanding on Table 3:
#           M   S   A
#   R1toR2  2   -   4
#   R1toR3  1   -   2
#   R2toR4  -   -   2

# R1toR2(X, Y, Z, Ta, Tb) = (X+Y, Y-X, Z+Z, (d*Ta*Tb) + (d*Ta*Tb))
def R1toR2(P):
    TaTb = GFp2.mul(P[3], P[4])
    dTaTb = GFp2.mul(d, TaTb)
    return (
        GFp2.add(P[0], P[1]),
        GFp2.sub(P[1], P[0]),
        GFp2.add(P[2], P[2]),
        GFp2.add(dTaTb, dTaTb)
    )

# R1toR3(X, Y, Z, Ta, Tb) = (X+Y, Y-X, Z, Ta*Tb)
def R1toR3(P):
    return(
        GFp2.add(P[0], P[1]),
        GFp2.sub(P[1], P[0]),
        P[2],
        GFp2.mul(P[3], P[4])
    )

# R2toR4(XpY, YmX, Z, T)  = (XpY-YmX, YmX-XmY, Z) // because projective
def R2toR4(P):
    return (
        GFp2.sub(P[0], P[1]),
        GFp2.add(P[1], P[0]),
        P[2]
    )

# Take a curve in any representation and clear the projective factor
def normalize(P):
    zi = GFp2.inv(P[2])
    xz = GFp2.mul(P[0], zi)
    yz = GFp2.mul(P[1], zi)
    if len(P) == 3: # R4
        return (xz, yz, one)
    if len(P) == 5: # R1
        taz = GFp2.mul(P[3], zi)
        tbz = GFp2.mul(P[4], zi)
        return (xz, yz, one, taz, tbz)
    raise Exception("Representation unsupported for normalization")

# Take a curve in affint form an put it in R1
def setup(P):
    return (P[0], P[1], one, P[0], P[1])

###

def reps_test():
    x  = (0, 1)
    x2 = (0, 2)
    y  = (2, 0)
    y2 = (4, 0)
    xy = (2, 1)
    yx = (2, p-1)
    z  = (3, 4)
    z2 = (6, 8)
    ta  = (5, 0)
    tb  = (1, 6)
    t   = (5, 30)
    td2 = GFp2.mul((2,0), GFp2.mul(d, t))

    r1 = (x, y, z, ta, tb)
    r2 = (xy, yx, z2, td2)
    r3 = (xy, yx, z, t)
    r4 = (x2, y2, z2)

    test("R1toR2", R1toR2(r1), r2)
    test("R1toR3", R1toR3(r1), r3)
    test("R2toR4", R2toR4(r2), r4)

    zi = GFp2.inv(z)
    xz = GFp2.mul(x, zi)
    yz = GFp2.mul(y, zi)
    taz = GFp2.mul(ta, zi)
    tbz = GFp2.mul(tb, zi)
    test("norm1", normalize(r1), (xz, yz, one, taz, tbz))
    test("norm4", normalize(r4), (xz, yz, one))

#reps_test()

########## Curve Arithmetic ##########

# ADD_core
#   R3 + R2 -> R1
#   7M    -S    4A
def ADD_core(P, Q):
    # Order of things is slightly cleaned up
    chi  = GFp2.mul(P[0], Q[0])
    xi   = GFp2.mul(P[1], Q[1])
    tau  = GFp2.mul(P[2], Q[2])
    zeta = GFp2.mul(P[3], Q[3])

    alpha = GFp2.add(tau, zeta)
    beta  = GFp2.sub(tau, zeta)
    mu    = GFp2.add(chi, xi)
    nu    = GFp2.sub(chi, xi)

    return (
        GFp2.mul(nu, beta),
        GFp2.mul(mu, alpha),
        GFp2.mul(alpha,  beta),
        mu,
        nu
    )

# ADD
#   R1 + R2 -> R1
#   8M    -S    6A
def ADD(P, Q):
    return ADD_core(R1toR3(P), Q)

# DBL
#   R4 -> R1
#   3M    4S    6A
def DBL(P):
    alpha = GFp2.sqr(P[0])        # alpha = X^2
    beta  = GFp2.sqr(P[1])        # beta = Y^2
    zeta  = GFp2.sqr(P[2])        # zeta = Z^2

    sigma = GFp2.add(P[0], P[1])  # sigma = X + Y
    tau   = GFp2.sqr(sigma)       # tau = sigma^2

    t1 = GFp2.add(alpha, beta)    # t1 = alpha + beta
    t2 = GFp2.sub(beta, alpha)    # t2 = beta - alpha
    t3 = GFp2.sub(tau, t1)        # t3 = tau  - t1
    t4 = GFp2.addsub(zeta, t2)    # t4 = zeta + zeta - t2

    return (
        GFp2.mul(t3, t4),         # X  = t3 * t4
        GFp2.mul(t1, t2),         # Y  = t1 * t2
        GFp2.mul(t2, t4),         # Z  = t2 * t4
        t3,                     # Ta = t3
        t1                      # Tb = t1
    )

###

def core_test():
    TEST_LOOPS = 1000

    # Test doubling
    A = setup(G)
    for i in range(TEST_LOOPS):
        A = DBL(A)
    A = normalize(A)
    doubleP = setup(((0x2C3FD8822C82270FC9099C54855859D6, 0x4DA5B9E83AA7A1B2A7B3F6E2043E8E68),
                     (0x2001EB3A576883963EE089F0EB49AA14, 0x0FFDB0D761421F501FEE5617A7E954CD)))
    testpt("double", A, doubleP)

    # Test that negation works
    P = setup(G)
    NP = setup((GFp2.neg(P[0]), P[1]))
    NP = R1toR2(NP)
    Z = normalize(ADD(P, NP))
    testpt("negation", Z, setup(O))

    # Test that the neutral element is neutral
    P = setup(G)
    Q = setup(O)
    PP = normalize(ADD(P, R1toR2(Q)))
    testpt("neutral-r", PP, P)
    P = setup(G)
    Q = setup(O)
    PP = normalize(ADD(Q, R1toR2(P)))
    testpt("neutral-l", PP, P)

    # Test point doubling by addition
    A = G
    P = setup(A)
    for i in range(TEST_LOOPS):
        Q = R1toR2(P)
        P = ADD(P, Q)
    A = normalize(P)
    testpt("double-add", A, doubleP)

    # Test repeated addition of the same point
    P = setup(G)
    Q = R1toR2(P)
    P = DBL(P)
    for i in range(TEST_LOOPS):
        P = ADD(P, Q)
    A = normalize(P)
    thouP = setup(((0x3E243958590C4D906480B1EF0A151DB0, 0x5327AF7D84238CD0AA270F644A65D473),
                   (0x3EF69A49CB7E02375E06003D73C43EB1, 0x293EB1E26DD23B4E4E752648AC2EF0AB)))
    testpt("addition", A, thouP)

#core_test()


########## Endomorphism ##########

## TODO: These are mechanical translations; clean them up

def tau(P):
    ctau1= (0x1964DE2C3AFAD20C74DCD57CEBCE74C3,
            0x000000000000000C0000000000000012)

    (Px, Py, Pz, Pta, Ptb) = P
    t0 = GFp2.sqr(Px)
    t1 = GFp2.sqr(Py)
    Px = GFp2.mul(Px, Py)
    Py = GFp2.sqr(Pz)
    Pz = GFp2.add(t0, t1)
    Py = GFp2.add(Py, Py)
    t0 = GFp2.sub(t0, t1)
    Py = GFp2.neg(Py)
    Px = GFp2.mul(Px, t0)
    Py = GFp2.sub(Py, t0)
    Px = GFp2.mul(Px, ctau1)
    Py = GFp2.mul(Py, Pz)
    Pz = GFp2.mul(Pz, t0)
    return (Px, Py, Pz, Pta, Ptb)


def tau_dual(P):
    ctaudual1 = (0x4AA740EB230586529ECAA6D9DECDF034,
                 0x7FFFFFFFFFFFFFF40000000000000011)

    (Px, Py, Pz, Pta, Ptb) = P
    t0  = GFp2.sqr(Px)
    Pta = GFp2.sqr(Pz)
    t1  = GFp2.sqr(Py)
    Pz  = GFp2.add(Pta, Pta)
    Pta = GFp2.sub(t1, t0)
    t0  = GFp2.add(t0, t1)
    Px  = GFp2.mul(Px, Py)
    Pz  = GFp2.sub(Pz, Pta)
    Ptb = GFp2.mul(Px, ctaudual1)
    Py  = GFp2.mul(Pz, Pta)
    Px  = GFp2.mul(Ptb, t0)
    Pz  = GFp2.mul(Pz, t0)
    return (Px, Py, Pz, Pta, Ptb)


def delphidel(P):
    cphi0 = (0x0000000000000005FFFFFFFFFFFFFFF7,
             0x2553A0759182C3294F65536CEF66F81A)
    cphi1 = (0x00000000000000050000000000000007,
             0x62C8CAA0C50C62CF334D90E9E28296F9)
    cphi2 = (0x000000000000000F0000000000000015,
             0x78DF262B6C9B5C982C2CB7154F1DF391)
    cphi3 = (0x00000000000000020000000000000003,
             0x5084C6491D76342A92440457A7962EA4)
    cphi4 = (0x00000000000000030000000000000003,
             0x12440457A7962EA4A1098C923AEC6855)
    cphi5 = (0x000000000000000A000000000000000F,
             0x459195418A18C59E669B21D3C5052DF3)
    cphi6 = (0x00000000000000120000000000000018,
             0x0B232A8314318B3CCD3643A78A0A5BE7)
    cphi7 = (0x00000000000000180000000000000023,
             0x3963BC1C99E2EA1A66C183035F48781A)
    cphi8 = (0x00000000000000AA00000000000000F0,
             0x1F529F860316CBE544E251582B5D0EF0)
    cphi9 = (0x00000000000008700000000000000BEF,
             0x0FD52E9CFE00375B014D3E48976E2505)

    (Px, Py, Pz, Pta, Ptb) = P
    t4 = GFp2.sqr(Pz)
    t3 = GFp2.mul(Py, Pz)
    t0 = GFp2.mul(t4, cphi4)
    t2 = GFp2.sqr(Py)
    t0 = GFp2.add(t0, t2)
    t1 = GFp2.mul(t3, cphi3)
    t5 = GFp2.sub(t0, t1)
    t0 = GFp2.add(t0, t1)
    t0 = GFp2.mul(t0, Pz)
    t1 = GFp2.mul(t3, cphi1)
    t0 = GFp2.mul(t0, t5)
    t5 = GFp2.mul(t4, cphi2)
    t5 = GFp2.add(t2, t5)
    t6 = GFp2.sub(t1, t5)
    t1 = GFp2.add(t1, t5)
    t6 = GFp2.mul(t6, t1)
    t6 = GFp2.mul(t6, cphi0)
    Px = GFp2.mul(Px, t6)
    t6 = GFp2.sqr(t2)
    t2 = GFp2.sqr(t3)
    t3 = GFp2.sqr(t4)
    t1 = GFp2.mul(t2, cphi8)
    t5 = GFp2.mul(t3, cphi9)
    t1 = GFp2.add(t1, t6)
    t2 = GFp2.mul(t2, cphi6)
    t3 = GFp2.mul(t3, cphi7)
    t1 = GFp2.add(t1, t5)
    t2 = GFp2.add(t2, t3)
    t1 = GFp2.mul(t1, Py)
    Py = GFp2.add(t6, t2)
    Px = GFp2.mul(Px, t1)
    Py = GFp2.mul(Py, cphi5)
    Px = GFp2.conj(Px)
    Py = GFp2.mul(Py, Pz)
    Pz = GFp2.mul(t0, t1)
    Py = GFp2.mul(Py, t0)
    Pz = GFp2.conj(Pz)
    Py = GFp2.conj(Py)
    return (Px, Py, Pz, Pta, Ptb)


def delpsidel(P):
    cpsi1 = (0x2AF99E9A83D54A02EDF07F4767E346EF,
             0x00000000000000DE000000000000013A)
    cpsi2 = (0x00000000000000E40000000000000143,
             0x21B8D07B99A81F034C7DEB770E03F372)
    cpsi3 = (0x00000000000000060000000000000009,
             0x4CB26F161D7D69063A6E6ABE75E73A61)
    cpsi4 = (0x7FFFFFFFFFFFFFF9FFFFFFFFFFFFFFF6,
             0x334D90E9E28296F9C59195418A18C59E)

    (Px, Py, Pz, Pta, Ptb) = P
    Px = GFp2.conj(Px)        # fpneg(Px[1]) # ?
    Py = GFp2.conj(Py)        # fpneg(Py[1]) # ?
    Pz = GFp2.conj(Pz)        # fpneg(Pz[1]) # ?
    t2 = GFp2.sqr(Pz)
    t0 = GFp2.sqr(Px)
    Px = GFp2.mul(Px, t2)
    Pz = GFp2.mul(t2, cpsi2)
    t1 = GFp2.mul(t2, cpsi3)
    t2 = GFp2.mul(t2, cpsi4)
    Pz = GFp2.add(t0, Pz)
    t2 = GFp2.add(t0, t2)
    t1 = GFp2.add(t0, t1)
    t2 = GFp2.neg(t2)
    Pz = GFp2.mul(Pz, Py)
    Px = GFp2.mul(Px, t2)
    Py = GFp2.mul(t1, Pz)
    Px = GFp2.mul(Px, cpsi1)
    Pz = GFp2.mul(Pz, t2)
    return (Px, Py, Pz, Pta, Ptb)

def phi(P):
    return tau_dual(delphidel(tau(P)))

def psi(P):
    return tau_dual(delpsidel(tau(P)))

def to64(x):
    return x & 0xffffffffffffffff

def decompose(k):
    ell1 = 0x0000000000000007FC5BB5C5EA2BE5DFF75682ACE6A6BD66259686E09D1A7D4FL
    ell2 = 0x00000000000000038FD4B04CAA6C0F8A2BD235580F468D8DD1BA1D84DD627AFBL
    ell3 = 0x0000000000000000D038BF8D0BFFBAF6C42BD6C965DCA9029B291A33678C203CL
    ell4 = 0x00000000000000031B073877A22D841081CBDC3714983D8212E5666B77E7FDC0L
    b11 = 0x0906FF27E0A0A196L
    b12 = 0x1363E862C22A2DA0L
    b13 = 0x07426031ECC8030FL
    b14 = 0x084F739986B9E651L
    b21 = 0x1D495BEA84FCC2D4L
    b24 = 0x25DBC5BC8DD167D0L
    b31 = 0x17ABAD1D231F0302L
    b32 = 0x02C4211AE388DA51L
    b33 = 0x2E4D21C98927C49FL
    b34 = 0x0A9E6F44C02ECD97L
    b41 = 0x136E340A9108C83FL
    b42 = 0x3122DF2DC3E0FF32L
    b43 = 0x068A49F02AA8A9B5L
    b44 = 0x18D5087896DE0AEAL
    c1  = 0x72482C5251A4559CL
    c2  = 0x59F95B0ADD276F6CL
    c3  = 0x7DD2D17C4625FA78L
    c4  = 0x6BC57DEF56CE8877L

    a1 = to64((k * ell1) >> 256)
    a2 = to64((k * ell2) >> 256)
    a3 = to64((k * ell3) >> 256)
    a4 = to64((k * ell4) >> 256)

    temp = to64(k - a1*b11 - a2*b21 - a3*b31  - a4*b41 + c1)
    mask = 0xffffffffffffffff * ((temp & 1) ^ 1)

    # XXX: This matches the code in FourQlib, but how does it line up with
    # the values in Proposition 5?
    return (
        to64(temp                                   + (mask & b41)),
        to64(a1*b12 + a2     - a3*b32 - a4*b42 + c2 + (mask & b42)),
        to64(a3*b33 - a1*b13 - a2     + a4*b43 + c3 - (mask & b43)),
        to64(a1*b14 - a2*b24 - a3*b34 + a4*b44 + c4 - (mask & b44))
    )

def recode(scalars):
    (a0, a1, a2, a3) = scalars
    sign_masks = range(65)
    digits = range(65)

    for i in range(64):
        a0 >>= 1
        bit0 = a0 & 1
        sign_masks[i] = CTRUE * (bit0)

        bit = a1 & 1
        carry = (bit0 | bit) ^ bit0
        a1 = (a1 >> 1) + carry
        digits[i] = bit

        bit = a2 & 1
        carry = (bit0 | bit) ^ bit0
        a2 = (a2 >> 1) + carry
        digits[i] += (bit << 1)

        bit = a3 & 1
        carry = (bit0 | bit) ^ bit0
        a3 = (a3 >> 1) + carry
        digits[i] += (bit << 2)

    sign_masks[64] = CTRUE
    digits[64] = a1 + (a2 << 1) + (a3 << 2)

    return digits, sign_masks

###

def endo_test():
    TEST_LOOPS = 1000

    # Test psi endomorphism
    P = setup(G)
    for i in range(TEST_LOOPS):
        P = psi(P)
    A = normalize(P)
    psiP = setup(((0x75AF54EDB41A2B93D8F3C8C24A2BC7E2, 0x065249F9EDE0C7984DE2466701F009A9),
                  (0x06DBB85BFFB7C21E1C6E119ADD608104, 0x060A30903424BF13FD234D6C4CFA3EC1)))
    testpt("psi", A, psiP)

    # Test phi endomorphism
    P = setup(G)
    for i in range(TEST_LOOPS):
        P = phi(P)
    A = normalize(P)
    phiP = setup(((0x5550AAB9E7A620EED5B5A3061287DB16, 0x3E61EBB9A1CB0210EC321E6CF33610FC),
                  (0x5474BF8EC55603AE7E2851D5A8E83FB9, 0x5476093DBF8BF6BFA5077613491788D5)))
    testpt("phi", A, phiP)

    # Test decomposition
    test_cases = [
        [0x92990788d66bf558052d112f5498111747b3e28c55984d43fed8c8822ad9f1a7L,
         (0xa8ea3f673f711e51L, 0xa08d1eae0b9e071dL, 0x55c8df690050276fL, 0x6396739dda88830f)],
        [0x48e5ca2a675ab49ca214b884813935024b0c61edc8d1305fe5230df341623348L,
         (0xa53ec4631945b875L, 0x521c0ba1261c1934L, 0x5c50ce912909185cL, 0x93b3c70960b44bad)],
        [0xae20e251c36cfa5be4d9f3d5a5edfed305a1e8f7f6394d9be58a15c4b0f1c5e9L,
         (0xa621ada9b3499c9fL, 0x7cd17e0095e7aae6L, 0x6e8d23b5bd10bb43L, 0x7f18c69f3025234c)],
        [0xb2c950abc87a55442cc00f1e3ac38f81b7e95036fd191ea134ff616d9806e10cL,
         (0x9b30a872ebea83afL, 0x8f6c73350447c9c3L, 0x72fdc76e3456d087L, 0x6ba39ba159b0c13d)],
        [0x8e2958a1475ed70762340e9797788e0061f21fcebd67889fdd4f4ce2b5f6b2deL,
         (0xbe8f3583a0934333L, 0xab45bf6d1bf80b37L, 0x4a19fc5cffe97809L, 0x5ea3baf1a1206442)],
    ]

    for t in test_cases:
        scalars = decompose(t[0])
        test("decompose", scalars, t[1])

    # Test recoding
    passed = True
    for n in range(TEST_LOOPS):
        k = getrandbits(256)
        scalars = decompose(k)
        digits, sign_masks = recode(scalars)

        a1 = a2 = a3 = a4 = 0
        for i in range(64, -1, -1):
            a1 *= 2
            a2 *= 2
            a3 *= 2
            a4 *= 2

            if sign_masks[i] > 0:
                a1 += 1
                a2 += (digits[i] >> 0) & 1
                a3 += (digits[i] >> 1) & 1
                a4 += (digits[i] >> 2) & 1
            else:
                a1 -= 1
                a2 -= (digits[i] >> 0) & 1
                a3 -= (digits[i] >> 1) & 1
                a4 -= (digits[i] >> 2) & 1

        if scalars[0] != a1 or scalars[1] != a2 or scalars[2] != a3 or scalars[3] != a4:
            passed = False
            break
    if passed:
        print "[PASS] recode"
    else:
        print "[FAIL] recode"

#endo_test()


########## Multiplication ##########

# Naive multiply, as a baseline
# R1 -> R1
def MUL_noendo(m, P):
    bits = [m & (1 << i) for i in range(256)]

    P2 = R1toR2(P)
    S = setup(O)

    for i in range(255, -1, -1):
        S = DBL(S)
        if bits[i] > 0:
            S = ADD(S, P2)
    return S

# See Algorithm 2
def MUL(m, P):
    # Compute endomorphisms
    Q = R1toR3(phi(P))
    R = R1toR3(psi(P))
    S = R1toR3(psi(phi(P)))

    # Precompute lookup table
    T = range(16)
    T[0] = R1toR2(P)
    T[1] = R1toR2(ADD_core(T[0], Q)) # T[1] = P+Q
    T[2] = R1toR2(ADD_core(T[0], R)) # T[2] = P+R
    T[3] = R1toR2(ADD_core(T[1], R)) # T[3] = P+Q+R
    T[4] = R1toR2(ADD_core(T[0], S)) # T[4] = P+S
    T[5] = R1toR2(ADD_core(T[4], Q)) # T[5] = P+Q+S
    T[6] = R1toR2(ADD_core(T[4], R)) # T[6] = P+R+S
    T[7] = R1toR2(ADD_core(T[6], Q)) # T[7] = P+Q+R+S

    # Scalar decomposition
    a = decompose(m)

    # Scalar recoding
    digits, sign_masks = recode(a)

    # Main loop
    def condneg(c, P):
        # Negates if c == CFALSE
        # In R2, neg(x, y, z, t) = (y, x, z, -t)
        return (
            GFp2.cswap(c, P[0], P[1]),
            GFp2.cswap(c, P[1], P[0]),
            P[2],
            GFp2.cswap(c, P[3], GFp2.neg(P[3]))
        )

    Q = condneg(sign_masks[64], T[digits[64]])
    Q = R2toR4(Q)
    for i in range(63, -1, -1):
        Q = DBL(Q)
        Q = ADD(Q, condneg(sign_masks[i], T[digits[i]]))
    return Q

def mul_test():
    TEST_LOOPS = 1000

    # Test multiplication by one and two (no endomorphisms)
    A = setup(G)
    B = MUL_noendo(1, A)
    testpt("mul-noendo-*1", B, A)
    A2 = DBL(A)
    B2 = MUL_noendo(2, A)
    testpt("mul-noendo-*2", B2, A2)

    # Test multiplication by one and two (with endomorphisms)
    A = setup(G)
    B = MUL(1, A)
    testpt("mul-*1", B, A)
    A2 = DBL(A)
    B2 = MUL(2, A)
    testpt("mul-*2", B2, A2)

    # Test multiply (no cofactor clearing, with/without endomorphism)
    # Also use this as a light performance test
    scalar = [0x3AD457AB55456230, 0x3A8B3C2C6FD86E0C, 0x7E38F7C9CFBB9166, 0x0028FD6CBDA458F0]
    coeff = range(TEST_LOOPS)
    for i in range(TEST_LOOPS):
         scalar[1] = scalar[2]
         scalar[2] += scalar[0]
         scalar[2] &= 0xffffffffffffffff
         coeff[i] = (scalar[0] << (0 * 64)) + (scalar[1] << (1 * 64)) + \
                    (scalar[2] << (2 * 64)) + (scalar[3] << (3 * 64))

    def testmul(label, mul):
        A = setup(G)
        for i in range(TEST_LOOPS):
            A = mul(coeff[i], A)

        mulP = setup(((0x257C122BBFC94A1BDFD2B477BD494BEF, 0x469BF80CB5B11F01769593547237C459),
                      (0x0901B3817C0E936C281C5067996F3344, 0x570B948EACACE2104FE8C429915F1245)))
        testpt(label, A, mulP)

    tic = time()
    testmul("mul-noendo", MUL_noendo)
    toc = time()
    elapsed_noendo = toc - tic

    tic = time()
    testmul("mul", MUL)
    toc = time()
    elapsed = toc - tic
    savings = 1.0 - elapsed * 1.0 / elapsed_noendo

    print "endo={:.2f} noendo={:.2f} savings={:2.1f}%".format(elapsed, elapsed_noendo, savings * 100)


mul_test()


########## Operation counts ##########

def opcount_test():
    opcounts = {}

    # R1toR2
    P = setup(G)
    GFp2.ctr_reset()
    Q = R1toR2(P)
    opcounts["R1toR2"] = GFp2.ctr()

    # R1toR3
    P = setup(G)
    GFp2.ctr_reset()
    Q = R1toR3(P)
    opcounts["R1toR3"] = GFp2.ctr()

    # R2toR4
    P = R1toR2(setup(G))
    GFp2.ctr_reset()
    Q = R2toR4(P)
    opcounts["R2toR4"] = GFp2.ctr()

    # ADD_core
    P = R1toR3(setup(G))
    Q = R1toR2(setup(G))
    GFp2.ctr_reset()
    R = ADD_core(P, Q)
    opcounts["ADD_core"] = GFp2.ctr()

    # ADD
    P = setup(G)
    Q = R1toR2(setup(G))
    GFp2.ctr_reset()
    R = ADD(P, Q)
    opcounts["ADD"] = GFp2.ctr()

    # DBL
    P = setup(G)
    GFp2.ctr_reset()
    Q = DBL(P)
    opcounts["DBL"] = GFp2.ctr()

    # MUL_noendo
    P = setup(G)
    m = p - 10
    GFp2.ctr_reset()
    Q = MUL_noendo(m, P)
    opcounts["MUL_noendo"] = GFp2.ctr()

    # Phi
    P = setup(G)
    GFp2.ctr_reset()
    phiP = phi(P)
    opcounts["phi"] = GFp2.ctr()

    # Psi
    P = setup(G)
    GFp2.ctr_reset()
    psiP = psi(P)
    opcounts["psi"] = GFp2.ctr()

    # Phi * Psi
    P = setup(G)
    GFp2.ctr_reset()
    psiphiP = psi(phiP)
    opcounts["psiphi"] = GFp2.ctr()

    # Preload
    opcounts["preload"] = (
        opcounts["phi"][0] + opcounts["psi"][0] + opcounts["psiphi"][0],
        opcounts["phi"][1] + opcounts["psi"][1] + opcounts["psiphi"][1],
        opcounts["phi"][2] + opcounts["psi"][2] + opcounts["psiphi"][2],
        opcounts["phi"][3] + opcounts["psi"][3] + opcounts["psiphi"][3]
    )

    # MUL
    m = getrandbits(256)
    P = setup(G)
    GFp2.ctr_reset()
    mP = MUL(m, P)
    opcounts["MUL"] = GFp2.ctr()

    rows = ["R1toR2", "R1toR3", "R2toR4", "ADD_core", "ADD", "DBL",
            "phi", "psi", "psiphi", "preload", "MUL", "MUL_noendo"]
    print "{:10s} {:>7s} {:>7s} {:>7s} {:>7s}".format("", "M", "S", "A", "I")
    for name in rows:
        if name not in opcounts:
            continue
        opctr = opcounts[name]
        print "{:10s} {:7.1f} {:7.1f} {:7.1f} {:7.1f}".format(name, opctr[2], opctr[1], opctr[0], opctr[3])

opcount_test()
