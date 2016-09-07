#!/usr/bin/env python

from fields import GFp25519 as GFp, p25519 as p

########## Constants ##########

# The curve parameter
d = 37095705934669439343138083508754565189542113879843219016388785533085940283555L

# The generator
Gx = 15112221349535400772501151409588531511454012693041857206046113283949847762202
Gy = 46316835694926478169428394003475163141307993866256225615783033603165251855960

########## Point encoding / decoding ##########

# BEGIN from RFC 7748
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

########## Point Multiplication ##########

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


def x25519(k, u):
    kn = decodeScalar(k)
    un = decodeUCoord(u)
    return encodeUCoord(x25519_inner(kn, un))

########## Self test for correctness ##########
########## (cases are from RFC 7748) ##########

def test_x25519():
    k0 = 'a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4'.decode('hex')
    u0 = 'e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c'.decode('hex')
    r0 = 'c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552'.decode('hex')

    rp = x25519(k0, u0)
    test.test('rfc-0', rp, r0)

    k1 = '4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d'.decode('hex')
    u1 = 'e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493'.decode('hex')
    r1 = '95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957'.decode('hex')
    test.test('rfc-1', x25519(k1, u1), r1)

    k = '0900000000000000000000000000000000000000000000000000000000000000'.decode('hex')
    u = '0900000000000000000000000000000000000000000000000000000000000000'.decode('hex')
    k1i = '422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079'.decode('hex')
    k1k = '684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51'.decode('hex')

    # XXX: Not currently tested
    k1m = '7c3911e0ab2586fd864497297e575e6f3bc601c0883c30df5f4dd2d24f665424'.decode('hex')

    TEST_LOOPS = 1001
    for i in range(TEST_LOOPS):
        r = x25519(k, u)
        u = k
        k = r

        if i == 0:
            test.test("rfc-iter-1", k, k1i)
        elif i == 999:
            test.test("rfc-iter-1k", k, k1k)
        elif i == 999999:
            # XXX: Not currently tested
            test.test("rfc-iter-1m", k, k1m)

def test_dh():
    a  = '77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a'.decode('hex')
    KA = '8520f0098930a754748b7ddcb43ef75a0dbf3a0d26381af4eba4a98eaa9b4e6a'.decode('hex')
    b  = '5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb'.decode('hex')
    KB = 'de9edb7d7b7dc1b4d35b61c2ece435373f8343c85b78674dadfc7e146f882b4f'.decode('hex')
    K  = '4a5d9d5ba4ce2de1728e3bf480350f25e07e21c947d19e3376f09b3c1e161742'.decode('hex')
    nine = '0900000000000000000000000000000000000000000000000000000000000000'.decode('hex')

    KAt = x25519(a, nine)
    KBt = x25519(b, nine)
    KABt = x25519(a, x25519(b, nine))
    KBAt = x25519(b, x25519(a, nine))

    test.test('DH-KA', KAt, KA)
    test.test('DH-KB', KBt, KB)
    test.test('DH-KAB', KABt, K)
    test.test('DH-KBA', KBAt, K)

if __name__ == "__main__":
    import test
    test_x25519()
    test_dh()
