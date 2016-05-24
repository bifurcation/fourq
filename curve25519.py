#!/usr/bin/env python

from time import time

# The characteristic of the field
p = (1 << 255) - 19

# The curve parameter
d = 37095705934669439343138083508754565189542113879843219016388785533085940283555L

# The generator
Gx = 15112221349535400772501151409588531511454012693041857206046113283949847762202
Gy = 46316835694926478169428394003475163141307993866256225615783033603165251855960

# A 512-bit mask for use in cswap
mask = (1 << 512) - 1

class GFp:
    ctr_enabled = True
    A = 0
    S = 0
    M = 0
    I = 0

    @staticmethod
    def ctr_reset():
        GFp.A = 0
        GFp.S = 0
        GFp.M = 0
        GFp.I = 0

    @staticmethod
    def ctr():
        return (GFp.A, GFp.S, GFp.M, GFp.I)

    @staticmethod
    def cswap(c, x, y):
        if GFp.ctr_enabled:
            # Counting each XOR as an add (-ish)
            GFp.A += 3
        temp = (mask * c) & (x ^ y)
        return (temp ^ x, temp ^ y)

    @staticmethod
    def add(x, y):
        if GFp.ctr_enabled:
            GFp.A += 1
        return (x + y) % p

    @staticmethod
    def sub(x, y):
        if GFp.ctr_enabled:
            GFp.A += 1
        return (x - y) % p

    @staticmethod
    def mul(x, y):
        if GFp.ctr_enabled:
            GFp.M += 1
        return (x * y) % p

    @staticmethod
    def sqr(x):
        if GFp.ctr_enabled:
            GFp.S += 1
        return (x * x) % p

    # Shamelessly copied from AGL's copy of djb's code
    # https://github.com/agl/curve25519-donna/blob/master/curve25519-donna.c#L774
    @staticmethod
    def inv(z):
        if GFp.ctr_enabled:
            GFp.I += 1
        orig_ctr_enabled = GFp.ctr_enabled
        GFp.ctr_enabled = False
        z2 = GFp.sqr(z)                   # 2
        t1 = GFp.sqr(z2)                  # 4
        t0 = GFp.sqr(t1)                  # 8
        z9 = GFp.mul(t0, z)               # 9
        z11 = GFp.mul(z9, z2)             # 11
        t0 = GFp.sqr(z11)                 # 22
        z2_5_0 = GFp.mul(t0, z9)          # 2^5 - 2^0 = 31 = 22 + 9

        t0 = GFp.sqr(z2_5_0)              # 2^6 - 2^1
        t1 = GFp.sqr(t0)                  # 2^7 - 2^2
        t0 = GFp.sqr(t1)                  # 2^8 - 2^3
        t1 = GFp.sqr(t0)                  # 2^9 - 2^4
        t0 = GFp.sqr(t1)                  # 2^10 - 2^5
        z2_10_0 = GFp.mul(t0, z2_5_0)     # 2^10 - 2^0

        t0 = GFp.sqr(z2_10_0)             # 2^11 - 2^1
        t1 = GFp.sqr(t0)                  # 2^12 - 2^2
        for i in range(2, 10, 2):       # 2^20 - 2^10
            t0 = GFp.sqr(t1)
            t1 = GFp.sqr(t0)
        z2_20_0 = GFp.mul(t1, z2_10_0)    # 2^20 - 2^0

        t0 = GFp.sqr(z2_20_0)             # 2^21 - 2^1
        t1 = GFp.sqr(t0)                  # 2^22 - 2^1
        for i in range(2, 20, 2):       # 2^40 - 2^20
            t0 = GFp.sqr(t1)
            t1 = GFp.sqr(t0)
        z2_40_0 = GFp.mul(t1, z2_20_0)    # 2^40 - 2^0

        t0 = GFp.sqr(z2_40_0)             # 2^41 - 2^1
        t1 = GFp.sqr(t0)                  # 2^42 - 2^2
        for i in range(2, 10, 2):       # 2^50 - 2^10
            t0 = GFp.sqr(t1)
            t1 = GFp.sqr(t0)
        z2_50_0 = GFp.mul(t1, z2_10_0)    # 2^50 - 2^0

        t0 = GFp.sqr(z2_50_0)             # 2^51 - 2^1
        t1 = GFp.sqr(t0)                  # 2^52 - 2^2
        for i in range(2, 50, 2):       # 2^100 - 2^50
            t0 = GFp.sqr(t1)
            t1 = GFp.sqr(t0)
        z2_100_0 = GFp.mul(t1, z2_50_0)   # 2^100 - 2^0

        t0 = GFp.sqr(z2_100_0)            # 2^101 - 2^1
        t1 = GFp.sqr(t0)                  # 2^102 - 2^2
        for i in range(2, 100, 2):      # 2^200 - 2^100
            t0 = GFp.sqr(t1)
            t1 = GFp.sqr(t0)
        z2_200_0 = GFp.mul(t1, z2_100_0)  # 2^200 - 2^0

        t0 = GFp.sqr(z2_200_0)            # 2^201 - 2^1
        t1 = GFp.sqr(t0)                  # 2^202 - 2^2
        for i in range(2, 50, 2):       # 2^250 - 2^50
            t0 = GFp.sqr(t1)
            t1 = GFp.sqr(t0)
        t0 = GFp.mul(t1, z2_50_0)         # 2^250 - 2^0

        t1 = GFp.sqr(t0)                  # 2^251 - 2^1
        t0 = GFp.sqr(t1)                  # 2^252 - 2^2
        t1 = GFp.sqr(t0)                  # 2^253 - 2^3
        t0 = GFp.sqr(t1)                  # 2^254 - 2^4
        t1 = GFp.sqr(t0)                  # 2^255 - 2^5
        t0 = GFp.mul(t1, z11)             # 2^255 - 21
        GFp.ctr_enabled = orig_ctr_enabled
        return t0

def transform(bits, a24, k, u):
    x1 = u
    x2 = 1
    z2 = 0
    x3 = u
    z3 = 1
    swap = 0

    for t in range(bits-1, -1, -1):
        kt = (k >> t) & 1
        swap ^= kt

        (x2, x3) = GFp.cswap(swap, x2, x3)
        (z2, z3) = GFp.cswap(swap, z2, z3)
        swap = kt

        A = GFp.add(x2, z2)
        AA = GFp.sqr(A)
        B = GFp.sub(x2, z2)
        BB = GFp.sqr(B)
        E = GFp.sub(AA, BB)

        C = GFp.add(x3, z3)
        D = GFp.sub(x3, z3)
        DA = GFp.mul(D, A)
        CB = GFp.mul(C, B)

        FF = GFp.sqr(GFp.add(DA, CB))
        GG = GFp.sqr(GFp.sub(DA, CB))

        x3 = FF
        z3 = GFp.mul(x1, GG)
        x2 = GFp.mul(AA, BB)
        z2 = GFp.mul(E, GFp.add(AA, GFp.mul(a24, E)))

    (x2, x3) = GFp.cswap(swap, x2, x3)
    (z2, z3) = GFp.cswap(swap, z2, z3)
    return GFp.mul(x2, GFp.inv(z2))

def x25519_inner(k, u):
    bits = 255
    a24 = 121665
    return transform(bits, a24, k, u)

# ----- BEGIN from RFC 7748 -----
def decodeLittleEndian(b, bits):
    return sum([b[i] << 8*i for i in range((bits+7)/8)])

def decodeScalar(k):
    k_list = [ord(b) for b in k]
    k_list[0] &= 248
    k_list[31] &= 127
    k_list[31] |= 64
    return decodeLittleEndian(k_list, 255)

def decodeUCoord(u):
    bits = 255
    u_list = [ord(b) for b in u]
    # Ignore any unused bits.
    if bits % 8:
        u_list[-1] &= (1<<(bits%8))-1
    return decodeLittleEndian(u_list, bits)

def encodeUCoord(u):
    bits = 255
    u = u % p
    return ''.join([chr((u >> 8*i) & 0xff)
                    for i in range((bits+7)/8)])
# ----- BEGIN from RFC 7748 -----

def x25519(k, u):
    kn = decodeScalar(k)
    un = decodeUCoord(u)
    return encodeUCoord(x25519_inner(kn, un))

###

def test(label, sample, ref):
    if sample == ref:
        print "[PASS] {}".format(label)
    else:
        print "[FAIL] {} {}".format(label, sample)

def x25519_test():
    k0 = 'a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4'.decode('hex')
    u0 = 'e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c'.decode('hex')
    r0 = 'c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552'.decode('hex')

    GFp.ctr_reset()
    rp = x25519(k0, u0)
    opctr = GFp.ctr()
    print "{:10s} {:>7s} {:>7s} {:>7s} {:>7s}".format("", "M", "S", "A", "I")
    print "{:10s} {:7.1f} {:7.1f} {:7.1f} {:7.1f}".format("x25519", opctr[2], opctr[1], opctr[0], opctr[3])
    print

    test('rfc-0', rp, r0)

    k1 = '4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d'.decode('hex')
    u1 = 'e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493'.decode('hex')
    r1 = '95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957'.decode('hex')
    test('rfc-1', x25519(k1, u1), r1)

    k = '0900000000000000000000000000000000000000000000000000000000000000'.decode('hex')
    u = '0900000000000000000000000000000000000000000000000000000000000000'.decode('hex')
    k1i = '422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079'.decode('hex')
    k1k = '684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51'.decode('hex')

    # XXX: Not currently tested
    k1m = '7c3911e0ab2586fd864497297e575e6f3bc601c0883c30df5f4dd2d24f665424'.decode('hex')

    tic = time()
    for i in range(1,1001):
        r = x25519(k, u)
        u = k
        k = r

        if i == 1:
            test("rfc-iter-1", k, k1i)
        elif i == 1000:
            test("rfc-iter-1k", k, k1k)
        elif i == 1000000:
            # XXX: Not currently tested
            test("rfc-iter-1m", k, k1m)
    toc = time()
    elapsed = toc - tic
    print "x25519={:.2f}".format(elapsed)

x25519_test()
