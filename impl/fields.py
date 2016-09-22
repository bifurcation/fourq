#!/usr/bin/env python

import struct

p1271 = (1 << 127) - 1      # Characteristic of 4Q base field
p25519 = (1 << 255) - 19    # Characteristic of x25519 field
mask = (1 << 512) - 1       # Mask for use in x25519 CSWAP

class GFp:
    ctr_enabled = True
    A = 0
    S = 0
    M = 0
    I = 0

    half = 0x40000000000000000000000000000000

    @staticmethod
    def ctr_reset():
        GFp.A = 0
        GFp.S = 0
        GFp.M = 0
        GFp.I = 0

    @staticmethod
    def ctr():
        return (GFp.A, GFp1271.S, GFp1271.M, GFp1271.I)

    @staticmethod
    def add(x, y):
        if GFp.ctr_enabled:
            GFp.A += 1
        return (x + y) % p1271

    @staticmethod
    def sub(x, y):
        if GFp.ctr_enabled:
            GFp.A += 1
        return (x - y) % p1271

    @staticmethod
    def mul(x, y):
        if GFp.ctr_enabled:
            GFp.M += 1
        return (x * y) % p1271

    @staticmethod
    def sqr(x):
        if GFp.ctr_enabled:
            GFp.S += 1
        return (x * x) % p1271

    @staticmethod
    def neg(x):
        if GFp.ctr_enabled:
            GFp.A += 1
        return (p1271 - x) % p1271

    @staticmethod
    def select(c, x, y):
        if GFp.ctr_enabled:
            # Counting each XOR as an add (-ish)
            GFp.A += 3
        return y ^ ((mask * c) & (x ^ y))

    @staticmethod
    def inv(x):
        t2 = GFp.sqr(x)           # 2
        t2 = GFp.mul(x, t2)       # 3
        t3 = GFp.sqr(t2)          # 6
        t3 = GFp.sqr(t3)          # 12
        t3 = GFp.mul(t2, t3)      # 15
        t4 = GFp.sqr(t3)          # 30
        t4 = GFp.sqr(t4)          # 60
        t4 = GFp.sqr(t4)          # 120
        t4 = GFp.sqr(t4)          # 240
        t4 = GFp.mul(t3, t4)      # 2^8 - 2^0
        t5 = GFp.sqr(t4)          # 2^9 - 2^1
        for i in range(7):        # 2^16 - 2^8
            t5 = GFp.sqr(t5)
        t5 = GFp.mul(t4, t5)      # 2^16 - 2^0
        t2 = GFp.sqr(t5)          # 2^17 - 2^1
        for i in range(15):       # 2^32 - 2^16
            t2 = GFp.sqr(t2)
        t2 = GFp.mul(t5, t2)      # 2^32 - 2^0
        t1 = GFp.sqr(t2)          # 2^33 - 2^1
        for i in range(31):       # 2^64 - 2^32
            t1 = GFp.sqr(t1)
        t1 = GFp.mul(t1, t2)      # 2^64 - 2^0
        for i in  range(32):      # 2^96 - 2^32
            t1 = GFp.sqr(t1)
        t1 = GFp.mul(t1, t2)      # 2^96 - 2^0
        for i in  range(16):      # 2^112 - 2^16
            t1 = GFp.sqr(t1)
        t1 = GFp.mul(t1, t5)      # 2^112 - 2^0
        for i in  range(8):       # 2^120 - 2^8
            t1 = GFp.sqr(t1)
        t1 = GFp.mul(t1, t4)      # 2^120 - 2^0
        for i in  range(4):       # 2^124 - 2^4
            t1 = GFp.sqr(t1)
        t1 = GFp.mul(t1, t3)      # 2^124 - 2^0
        t1 = GFp.sqr(t1)          # 2^125 - 2^1
        t1 = GFp.mul(t1, x)       # 2^125 - 2^0
        t1 = GFp.sqr(t1)          # 2^126 - 2^1
        t1 = GFp.sqr(t1)          # 2^127 - 2^2
        return GFp.mul(t1, x)     # 2^127 - 3

    # 1/sqrt(x) = x^((p-3)/4) = x^(2^125 - 1)
    @staticmethod
    def invsqrt(x):
        xp = GFp.mul(x, x)                  # 2
        xp = GFp.mul(xp, xp)                # 4
        xp = GFp.mul(xp, x)                 # 5
        xp = GFp.mul(GFp.mul(xp, xp), xp)   # 15
        xp = GFp.mul(GFp.mul(xp, xp), x)    # 31 = 2^5 - 1

        accum = xp
        for i in range(24):
            for j in range(5):
                xp = GFp.sqr(xp)            # 2^(5(i+1)) - 2^(5i)
            accum = GFp.mul(xp, accum)      # 2^(5(i+1)) - 1
        return accum

    @staticmethod
    def toLittleEndian(x):
        return bytearray(struct.pack('<QQ', x % (1 << 64), x >> 64))

    @staticmethod
    def fromLittleEndian(x):
        x[15] &= 0x7F
        (lo, hi) = struct.unpack('<QQ', x)
        return (hi << 64) + lo

class GFp2:
    A = 0
    S = 0
    M = 0
    I = 0

    zero = (0, 0)
    one  = (1, 0)
    two  = (2, 0)
    half = (GFp.half, 0)

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
        return ((a[0] + b[0]) % p1271, (a[1] + b[1]) % p1271)

    @staticmethod
    def sub(a, b):
        GFp2.A += 1
        return ((a[0] - b[0]) % p1271, (a[1] - b[1]) % p1271)

    @staticmethod
    def mul(a, b):
        GFp2.M += 1
        a0b0 = (a[0] * b[0])
        a1b0 = (a[1] * b[0])
        a0b1 = (a[0] * b[1])
        a1b1 = (a[1] * b[1])
        return ((a0b0 - a1b1) % p1271, (a0b1 + a1b0) % p1271)

    @staticmethod
    def sqr(a):
        GFp2.S += 1
        a02   = (a[0] * a[0])
        a12   = (a[1] * a[1])
        a0a12 = (2 * a[0] * a[1])
        return ((a02 - a12) % p1271, a0a12 % p1271)

    @staticmethod
    def neg(a):
        GFp2.A += 1
        return ((p1271 - a[0]) % p1271, (p1271 - a[1]) % p1271)

    @staticmethod
    def conj(a):
        GFp2.A += 0.5
        return (a[0], (p1271 - a[1]) % p1271)

    @staticmethod
    def inv(a):
        GFp2.I += 1
        invmag = GFp.inv(GFp.add(GFp.sqr(a[0]), GFp.sqr(a[1])))
        GFp2.M -= 1 # To avoid double-counting
        GFp2.A -= .5
        return GFp2.mul((invmag, 0), GFp2.conj(a))

    # NB: Not constant-time
    @staticmethod
    def invsqrt(a):
        if a[1] == 0:
            t = GFp.invsqrt(a[0])
            if GFp.mul(a[0], GFp.sqr(t)) == 1:
                return (t, 0)
            else:
                return (0, t)

        def chi(x):
            return GFp.sqr(GFp.sqrt(x))

        n = GFp.add(GFp.sqr(a[0]), GFp.sqr(a[1]))
        s = GFp.invsqrt(n)
        c = GFp.mul(n, s)
        if GFp.mul(c, s) == -1:
            raise Exception('not square')

        delta = GFp.mul(GFp.add(a[0], c), GFp.half)
        g = GFp.invsqrt(delta)
        h = GFp.mul(delta, g)
        if GFp.mul(h, g) == -1:
            delta = GFp.mul(GFp.sub(a[0], c), GFp.half)
            g = GFp.invsqrt(delta)
            h = GFp.mul(delta, g)

        x0 = GFp.mul(h, s)
        x1 = GFp.neg(GFp.mul(GFp.mul(GFp.mul(a[1], s), g), GFp.half))
        return (x0, x1)

    # switch (c) {
    #   case 1:  return x;
    #   case 0: return y;
    # }
    @staticmethod
    def select(c, x, y):
        return (GFp.select(c, x[0], y[0]), GFp.select(c, x[1], y[1]))

class GFp25519:
    ctr_enabled = True
    A = 0
    S = 0
    M = 0
    I = 0

    @staticmethod
    def ctr_reset():
        GFp25519.A = 0
        GFp25519.S = 0
        GFp25519.M = 0
        GFp25519.I = 0

    @staticmethod
    def ctr():
        return (GFp25519.A, GFp25519.S, GFp25519.M, GFp25519.I)

    @staticmethod
    def cswap(c, x, y):
        if GFp25519.ctr_enabled:
            # Counting each XOR as an add (-ish)
            GFp25519.A += 3
        temp = (mask * c) & (x ^ y)
        return (temp ^ x, temp ^ y)

    @staticmethod
    def add(x, y):
        if GFp25519.ctr_enabled:
            GFp25519.A += 1
        return (x + y) % p25519

    @staticmethod
    def sub(x, y):
        if GFp25519.ctr_enabled:
            GFp25519.A += 1
        return (x - y) % p25519

    @staticmethod
    def mul(x, y):
        if GFp25519.ctr_enabled:
            GFp25519.M += 1
        return (x * y) % p25519

    @staticmethod
    def sqr(x):
        if GFp25519.ctr_enabled:
            GFp25519.S += 1
        return (x * x) % p25519

    # Copied from AGL's copy of djb's code
    # https://github.com/agl/curve25519-donna/blob/master/curve25519-donna.c#L774
    @staticmethod
    def inv(z):
        if GFp25519.ctr_enabled:
            GFp25519.I += 1
        orig_ctr_enabled = GFp25519.ctr_enabled
        GFp25519.ctr_enabled = False
        z2 = GFp25519.sqr(z)                   # 2
        t1 = GFp25519.sqr(z2)                  # 4
        t0 = GFp25519.sqr(t1)                  # 8
        z9 = GFp25519.mul(t0, z)               # 9
        z11 = GFp25519.mul(z9, z2)             # 11
        t0 = GFp25519.sqr(z11)                 # 22
        z2_5_0 = GFp25519.mul(t0, z9)          # 2^5 - 2^0 = 31 = 22 + 9

        t0 = GFp25519.sqr(z2_5_0)              # 2^6 - 2^1
        t1 = GFp25519.sqr(t0)                  # 2^7 - 2^2
        t0 = GFp25519.sqr(t1)                  # 2^8 - 2^3
        t1 = GFp25519.sqr(t0)                  # 2^9 - 2^4
        t0 = GFp25519.sqr(t1)                  # 2^10 - 2^5
        z2_10_0 = GFp25519.mul(t0, z2_5_0)     # 2^10 - 2^0

        t0 = GFp25519.sqr(z2_10_0)             # 2^11 - 2^1
        t1 = GFp25519.sqr(t0)                  # 2^12 - 2^2
        for i in range(2, 10, 2):              # 2^20 - 2^10
            t0 = GFp25519.sqr(t1)
            t1 = GFp25519.sqr(t0)
        z2_20_0 = GFp25519.mul(t1, z2_10_0)    # 2^20 - 2^0

        t0 = GFp25519.sqr(z2_20_0)             # 2^21 - 2^1
        t1 = GFp25519.sqr(t0)                  # 2^22 - 2^1
        for i in range(2, 20, 2):              # 2^40 - 2^20
            t0 = GFp25519.sqr(t1)
            t1 = GFp25519.sqr(t0)
        z2_40_0 = GFp25519.mul(t1, z2_20_0)    # 2^40 - 2^0

        t0 = GFp25519.sqr(z2_40_0)             # 2^41 - 2^1
        t1 = GFp25519.sqr(t0)                  # 2^42 - 2^2
        for i in range(2, 10, 2):              # 2^50 - 2^10
            t0 = GFp25519.sqr(t1)
            t1 = GFp25519.sqr(t0)
        z2_50_0 = GFp25519.mul(t1, z2_10_0)    # 2^50 - 2^0

        t0 = GFp25519.sqr(z2_50_0)             # 2^51 - 2^1
        t1 = GFp25519.sqr(t0)                  # 2^52 - 2^2
        for i in range(2, 50, 2):              # 2^100 - 2^50
            t0 = GFp25519.sqr(t1)
            t1 = GFp25519.sqr(t0)
        z2_100_0 = GFp25519.mul(t1, z2_50_0)   # 2^100 - 2^0

        t0 = GFp25519.sqr(z2_100_0)            # 2^101 - 2^1
        t1 = GFp25519.sqr(t0)                  # 2^102 - 2^2
        for i in range(2, 100, 2):             # 2^200 - 2^100
            t0 = GFp25519.sqr(t1)
            t1 = GFp25519.sqr(t0)
        z2_200_0 = GFp25519.mul(t1, z2_100_0)  # 2^200 - 2^0

        t0 = GFp25519.sqr(z2_200_0)            # 2^201 - 2^1
        t1 = GFp25519.sqr(t0)                  # 2^202 - 2^2
        for i in range(2, 50, 2):              # 2^250 - 2^50
            t0 = GFp25519.sqr(t1)
            t1 = GFp25519.sqr(t0)
        t0 = GFp25519.mul(t1, z2_50_0)         # 2^250 - 2^0

        t1 = GFp25519.sqr(t0)                  # 2^251 - 2^1
        t0 = GFp25519.sqr(t1)                  # 2^252 - 2^2
        t1 = GFp25519.sqr(t0)                  # 2^253 - 2^3
        t0 = GFp25519.sqr(t1)                  # 2^254 - 2^4
        t1 = GFp25519.sqr(t0)                  # 2^255 - 2^5
        t0 = GFp25519.mul(t1, z11)             # 2^255 - 21
        GFp25519.ctr_enabled = orig_ctr_enabled
        return t0

########## Self test for correctness ##########

def test_GFp():
    inv13 = GFp.inv(13)
    test.test("inv-1271", GFp.mul(inv13, 13), 1)

    invsqrt13 = GFp.invsqrt(13)
    test.test("invsqrt-1271", GFp.mul(13, GFp.sqr(invsqrt13)), 1)

def test_GFp2():
    one = (1,0)
    i = (0,1)
    x23 = (2, 3)
    x57 = (5, 7)

    test.test("1+i", GFp2.add(one, i), (1, 1))
    test.test("1*i", GFp2.mul(one, i), (0, 1))
    test.test("i*i", GFp2.mul(i, i), (p1271 - 1, 0))

    test.test("add",     GFp2.add(x23, x57), (7, 10))
    test.test("sub-pos", GFp2.sub(x57, x23), (3, 4))
    test.test("sub-neg", GFp2.sub(x23, x57), (p1271 - 3, p1271 - 4))
    test.test("mul",     GFp2.mul(x23, x57), (p1271 - 11, 29))
    test.test("sqr",     GFp2.sqr(x23), (p1271 - 5, 12))

    x23c = GFp2.conj(x23)
    x23i = GFp2.inv(x23)
    x23is = GFp2.invsqrt(x23)
    test.test("conj", x23c, (2, GFp.neg(3)))
    test.test("inv-1271-2", GFp2.mul(x23, x23i), one)
    test.test("invsqrt-1271-2", GFp2.mul(x23, GFp2.sqr(x23is)), one)

    select1 = GFp2.select(1, x23, x57)
    select0 = GFp2.select(0, x23, x57)
    test.test("select(1)", select1, x23)
    test.test("select(0)", select0, x57)

def test_GFp25519():
    inv13 = GFp25519.inv(13)
    test.test("inv-25519", GFp25519.mul(inv13, 13), 1)

if __name__ == "__main__":
    import test
    test_GFp()
    test_GFp2()
    test_GFp25519()
