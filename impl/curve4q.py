#!/usr/bin/env python

from random import getrandbits
from fields import GFp, GFp2, p1271

########## Definitions ##########

# Curve parameter as a field element tuple and its double
d = (0xe40000000000000142, 0x5e472f846657e0fcb3821488f1fc0c8d)

# Order of the curve group
N = 0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe0f79992fb2540ec7768ce7

# Neutral point coordinates
Ox = (0, 0)
Oy = (1, 0)

# Base point coordinates
Gx = (0x1A3472237C2FB305286592AD7B3833AA, 0x1E1F553F2878AA9C96869FB360AC77F6)
Gy = (0x0E3FEE9BA120785AB924A2462BCBB287, 0x6E1C4AF8630E024249A7C344844C8B5C)

# E: -x^2 + y^2 = 1 + d * x^2 * y^2
def PointOnCurve(P):
    (X, Y) = P
    X2 = GFp2.sqr(X)
    Y2 = GFp2.sqr(Y)
    LHS = GFp2.sub(Y2, X2)
    RHS = GFp2.add(GFp2.one, GFp2.mul(GFp2.mul(d, X2), Y2))
    return LHS == RHS

########## Point encoding / decoding ##########

def sign(X):
    s0 = X[0] >> 126
    s1 = X[1] >> 126
    if X[0] != 0:
        return s0
    else:
        return s1

def encode(X, Y):
    y0 = GFp.toLittleEndian(Y[0])
    y1 = GFp.toLittleEndian(Y[1])
    s = sign(X)
    y1[15] |= (s << 7)
    return y0 + y1

# NB: Not constant-time
def decode(B):
    if len(B) != 32:
        raise Exception("Malformed point: length {} != 32".format(len(B)))
    if B[15] & 0x80 != 0x00:
        raise Exception("Malformed point: reserved bit is not zero")

    s = B[31] >> 7
    B[31] &= 0x7F

    y0 = GFp.fromLittleEndian(B[:16])
    y1 = GFp.fromLittleEndian(B[16:])

    if y0 >= p1271 or y1 >= p1271:
        raise Exception("Malformed point: reserved bit is not zero")

    y = (y0, y1)
    y2 = GFp2.sqr(y)
    (u0, u1) = GFp2.sub(y2, GFp2.one)
    (v0, v1) = GFp2.add(GFp2.mul(d, y2), GFp2.one)

    t0 = GFp.add(GFp.mul(u0, v0), GFp.mul(u1, v1))
    t1 = GFp.sub(GFp.mul(u1, v0), GFp.mul(u0, v1))
    t2 = GFp.add(GFp.sqr(v0), GFp.sqr(v1))
    t3 = GFp.add(GFp.sqr(t0), GFp.sqr(t1))
    t3 = GFp.mul(GFp.invsqrt(t3), t3)

    t = GFp.mul(2, GFp.add(t0, t3))
    if t == 0:
        t = GFp.mul(GFp.two, GFp.sub(t0, t3))

    a = GFp.invsqrt(GFp.mul(t, GFp.mul(t2, GFp.sqr(t2))))
    b = GFp.mul(GFp.mul(a, t2), t)

    x0 = GFp.mul(b, GFp.half)
    x1 = GFp.mul(GFp.mul(a, t2), t1)
    if t != GFp.mul(t2, GFp.sqr(b)):
        x0, x1 = x1, x0

    x = (x0, x1)
    if sign(x) != s:
        x = GFp2.neg(x)

    if not PointOnCurve((x, y)):
        x = GFp2.conj(x)
    if not PointOnCurve((x, y)):
        raise Exception("Point not on curve")

    return (x, y)

########## Alternative Point Representations and Addition Laws ##########

def AffineToR1(X, Y):
    return (X, Y, GFp2.one, X, Y)

def R1toAffine(P):
    (X, Y, Z, Ta, Tb) = P
    Zi = GFp2.inv(Z)
    return (GFp2.mul(X, Zi), GFp2.mul(Y, Zi))

# R1toR2(X, Y, Z, Ta, Tb) = (X+Y, Y-X, Z+Z, 2d*Ta*Tb)
def R1toR2(P):
    (X, Y, Z, Ta, Tb) = P
    return (
        GFp2.add(X, Y),
        GFp2.sub(Y, X),
        GFp2.add(Z, Z),
        GFp2.mul(GFp2.mul(GFp2.two, d), GFp2.mul(Ta, Tb))
    )

# R1toR3(X, Y, Z, Ta, Tb) = (X+Y, Y-X, Z, Ta*Tb)
def R1toR3(P):
    (X, Y, Z, Ta, Tb) = P
    return(
        GFp2.add(X, Y),
        GFp2.sub(Y, X),
        Z,
        GFp2.mul(Ta, Tb)
    )

# R2toR4(N, D, E, F)  = (N-D, D-N, Z)
def R2toR4(P):
    (N, D, E, F) = P
    return (
        GFp2.sub(N, D),
        GFp2.add(D, N),
        E
    )

# R1/R4 -> R1
def DBL(P):
    (X1, Y1, Z1) = P[:3]
    A = GFp2.sqr(X1)
    B = GFp2.sqr(Y1)
    C = GFp2.mul(GFp2.two, GFp2.sqr(Z1))
    D = GFp2.add(A, B)
    E = GFp2.sub(GFp2.sqr(GFp2.add(X1, Y1)), D)
    F = GFp2.sub(B, A)
    G = GFp2.sub(C, F)
    X3 = GFp2.mul(E, G)
    Y3 = GFp2.mul(D, F)
    Z3 = GFp2. mul(F, G)
    Ta3 = E
    Tb3 = D
    return (X3, Y3, Z3, Ta3, Tb3)

# R3 + R2 -> R1
def ADD_core(P, Q):
    (N1, D1, E1, F1) = P
    (N2, D2, Z2, T2) = Q
    A = GFp2.mul(D1, D2)
    B = GFp2.mul(N1, N2)
    C = GFp2.mul(T2, F1)
    D = GFp2.mul(Z2, E1)
    E = GFp2.sub(B, A)
    F = GFp2.sub(D, C)
    G = GFp2.add(D, C)
    H = GFp2.add(B, A)
    X3 = GFp2.mul(E, F)
    Y3 = GFp2.mul(G, H)
    Z3 = GFp2.mul(F, G)
    Ta3 = E
    Tb3 = H
    return (X3, Y3, Z3, Ta3, Tb3)

# R1 + R2 -> R1
def ADD(P, Q):
    return ADD_core(R1toR3(P), Q)

########## Multiplication without Endomorphisms ##########

def table_windowed(P):
    Q = DBL(P)
    T = range(8)
    T[0] = R1toR2(P)
    for i in range(1, 8):
        T[i] = R1toR2(ADD(Q, T[i-1]))
    return T

# R1 -> R1
def MUL_windowed(m, P, table=None):
    # Check that P is in R1
    (X, Y, Z, Ta, Tb) = P

    # Negation in R2
    def R2neg(P):
        (N, D, E, F) = P
        return (D, N, E, GFp2.neg(F))

    # Constant-time point selection (in R2)
    def selectpt(c, P1, P2):
        (N1, D1, E1, F1) = P1
        (N2, D2, E2, F2) = P2
        return (
            GFp2.select(c, N1, N2),
            GFp2.select(c, D1, D2),
            GFp2.select(c, E1, E2),
            GFp2.select(c, F1, F2),
        )

    # Pre-compute [1]P, [3]P, ..., [15]P and negatives
    # XXX: Should check whether this table is for this point
    T = table
    if not T:
        T = table_windowed(P)
    nT = [R2neg(P) for P in T]

    # Pre-compute scalars
    d = range(63)
    reduced = m % N
    if reduced % 2 == 0:
        reduced += N
    for i in range(63):
        d[i] = (reduced % 32) - 16
        reduced = (reduced - d[i]) / 16
    d[62] = reduced
    ind = [(abs(di) - 1) / 2 for di in d]
    sgn = [di / abs(di) for di in d]
    sgn = [(s + 1) / 2 for s in sgn] # -1/+1 to 0/1

    # Compute the multiplication
    Q = R2toR4(selectpt(sgn[62], T[ind[62]], nT[ind[62]]))
    for i in range(61, -1, -1):
        Q = DBL(DBL(DBL(DBL(Q))))
        S = selectpt(sgn[i], T[ind[i]], nT[ind[i]])
        Q = ADD(Q, S)

    return Q


########## Endomorphisms ##########

ctau     = (0x1964de2c3afad20c74dcd57cebce74c3, 0x000000000000000c0000000000000012)
ctaudual = (0x4aa740eb230586529ecaa6d9decdf034, 0x7ffffffffffffff40000000000000011)

cphi0 = (0x0000000000000005fffffffffffffff7, 0x2553a0759182c3294f65536cef66f81a)
cphi1 = (0x00000000000000050000000000000007, 0x62c8caa0c50c62cf334d90e9e28296f9)
cphi2 = (0x000000000000000f0000000000000015, 0x78df262b6c9b5c982c2cb7154f1df391)
cphi3 = (0x00000000000000020000000000000003, 0x5084c6491d76342a92440457a7962ea4)
cphi4 = (0x00000000000000030000000000000003, 0x12440457a7962ea4a1098c923aec6855)
cphi5 = (0x000000000000000a000000000000000f, 0x459195418a18c59e669b21d3c5052df3)
cphi6 = (0x00000000000000120000000000000018, 0x0b232a8314318b3ccd3643a78a0a5be7)
cphi7 = (0x00000000000000180000000000000023, 0x3963bc1c99e2ea1a66c183035f48781a)
cphi8 = (0x00000000000000aa00000000000000f0, 0x1f529f860316cbe544e251582b5d0ef0)
cphi9 = (0x00000000000008700000000000000bef, 0x0fd52e9cfe00375b014d3e48976e2505)
cpsi1 = (0x2af99e9a83d54a02edf07f4767e346ef, 0x00000000000000de000000000000013a)
cpsi2 = (0x00000000000000e40000000000000143, 0x21b8d07b99a81f034c7deb770e03f372)
cpsi3 = (0x00000000000000060000000000000009, 0x4cb26f161d7d69063a6e6abe75e73a61)
cpsi4 = (0x7ffffffffffffff9fffffffffffffff6, 0x334d90e9e28296f9c59195418a18c59e)

def tau(P):
    (X1, Y1, Z1) = P
    A = GFp2.sqr(X1)
    B = GFp2.sqr(Y1)
    C = GFp2.add(A, B)
    D = GFp2.sub(A, B)
    X2 = GFp2.mul(GFp2.mul(GFp2.mul(ctau, X1), Y1), D)
    Y2 = GFp2.neg(GFp2.mul(GFp2.add(GFp2.mul(GFp2.two, GFp2.sqr(Z1)), D), C))
    Z2 = GFp2.mul(C, D)
    return (X2, Y2, Z2)

def tau_dual(P):
    (X1, Y1, Z1) = P
    A = GFp2.sqr(X1)
    B = GFp2.sqr(Y1)
    C = GFp2.add(A, B)
    Ta2 = GFp2.sub(B, A)
    D = GFp2.sub(GFp2.mul(GFp2.two, GFp2.sqr(Z1)), Ta2)
    Tb2 = GFp2.mul(GFp2.mul(ctaudual, X1), Y1)
    X2 = GFp2.mul(Tb2, C)
    Y2 = GFp2.mul(Ta2, D)
    Z2 = GFp2.mul(C, D)
    return (X2, Y2, Z2, Ta2, Tb2)

def upsilon(P):
    (X1, Y1, Z1) = P
    A = GFp2.mul(GFp2.mul(cphi0, X1), Y1)
    B = GFp2.mul(Y1, Z1)
    C = GFp2.sqr(Y1)
    D = GFp2.sqr(Z1)
    F = GFp2.sqr(D)
    G = GFp2.sqr(B)
    H = GFp2.sqr(C)
    I = GFp2.mul(cphi1, B)
    J = GFp2.add(C, GFp2.mul(cphi2, D))
    K = GFp2.add(GFp2.add(GFp2.mul(cphi8, G), H), GFp2.mul(cphi9, F))
    X2 = GFp2.mul(GFp2.add(I, J), GFp2.sub(I, J))
    X2 = GFp2.conj(GFp2.mul(GFp2.mul(A, K), X2))
    L = GFp2.add(C, GFp2.mul(cphi4, D))
    M = GFp2.mul(cphi3, B)
    N = GFp2.mul(GFp2.add(L, M), GFp2.sub(L, M))
    Y2 = GFp2.add(GFp2.add(H, GFp2.mul(cphi6, G)), GFp2.mul(cphi7, F))
    Y2 = GFp2.conj(GFp2.mul(GFp2.mul(GFp2.mul(cphi5, D), N), Y2))
    Z2 = GFp2.conj(GFp2.mul(GFp2.mul(B, K), N))
    return (X2, Y2, Z2)

def chi(P):
    (X1, Y1, Z1) = P
    A = GFp2.conj(X1)
    B = GFp2.conj(Y1)
    C = GFp2.sqr(GFp2.conj(Z1))
    D = GFp2.sqr(A)
    F = GFp2.sqr(B)
    G = GFp2.mul(B, GFp2.add(D, GFp2.mul(cpsi2, C)))
    H = GFp2.neg(GFp2.add(D, GFp2.mul(cpsi4, C)))
    X2 = GFp2.mul(GFp2.mul(GFp2.mul(cpsi1, A), C), H)
    Y2 = GFp2.mul(G, GFp2.add(D, GFp2.mul(cpsi3, C)))
    Z2 = GFp2.mul(G, H)
    return (X2, Y2, Z2)

def phi(P):
    return tau_dual(upsilon(tau(P[:3])))

def psi(P):
    return tau_dual(chi(tau(P[:3])))

########## Recoding ##########

b1 = [0x0906ff27e0a0a196, -0x1363e862c22a2da0,  0x07426031ecc8030f, -0x084f739986b9e651]
b2 = [0x1d495bea84fcc2d4, -0x0000000000000001,  0x0000000000000001,  0x25dbc5bc8dd167d0]
b3 = [0x17abad1d231f0302,  0x02c4211ae388da51, -0x2e4d21c98927c49f,  0x0a9e6f44c02ecd97]
b4 = [0x136e340a9108c83f,  0x3122df2dc3e0ff32, -0x068a49f02aa8a9b5, -0x18d5087896de0aea]

L1 = 0x7fc5bb5c5ea2be5dff75682ace6a6bd66259686e09d1a7d4f
L2 = 0x38fd4b04caa6c0f8a2bd235580f468d8dd1ba1d84dd627afb
L3 = 0x0d038bf8d0bffbaf6c42bd6c965dca9029b291a33678c203c
L4 = 0x31b073877a22d841081cbdc3714983d8212e5666b77e7fdc0

c  = [5 * b2[i] - 3 * b3[i] + 2 * b4[i] for i in range(4)]
cp = [c[i] + b4[i] for i in range(4)]

def decompose(m):
    def select(s, x, y):
        mask = (1 << 64) - 1
        return y ^ ((mask * s) & (x ^ y))

    t1 = (L1 * m) >> 256
    t2 = (L2 * m) >> 256
    t3 = (L3 * m) >> 256
    t4 = (L4 * m) >> 256

    a   = [m, 0, 0, 0]
    a   = [a[i] - t1*b1[i] - t2*b2[i] - t3*b3[i] - t4*b4[i] for i in range(4)]
    ac  = [a[i] + c[i]  for i in range(4)]
    acp = [a[i] + cp[i] for i in range(4)]

    s = ac[0] % 2
    v = [select(s, ac[i], acp[i]) for i in range(4)]
    return v

def recode(v):
    (v1, v2, v3, v4) = v
    vv = [v1, v2, v3, v4]

    def bit(x, n):
        return (x >> n) & 1

    d = range(65)
    m = range(65)
    for i in range(64):
        b1 = bit(vv[0], i+1)
        d[i] = 0
        m[i] = b1

        for j in [1, 2, 3]:
            bj = bit(vv[j], 0)
            d[i] += bj << (j - 1)
            c = (b1 | bj) ^ b1
            vv[j] = (vv[j] >> 1) + c

    d[64] = vv[1] + 2*vv[2] + 4*vv[3]
    m[64] = 1
    return (m, d)


########## Optimized multiplication ##########

def table_endo(P):
    Q = phi(P)
    R = psi(P)
    S = psi(Q)

    Q = R1toR3(Q)
    R = R1toR3(R)
    S = R1toR3(S)

    T = range(8)
    T[0] = R1toR2(P)                 # P
    T[1] = R1toR2(ADD_core(Q, T[0])) # P + Q
    T[2] = R1toR2(ADD_core(R, T[0])) # P + R
    T[3] = R1toR2(ADD_core(R, T[1])) # P + Q + R
    T[4] = R1toR2(ADD_core(S, T[0])) # P + S
    T[5] = R1toR2(ADD_core(S, T[1])) # P + Q + S
    T[6] = R1toR2(ADD_core(S, T[2])) # P + R + S
    T[7] = R1toR2(ADD_core(S, T[3])) # P + Q + R + S
    return T

def MUL_endo(m, P, table=None):
    # Check that P is in R1
    (X, Y, Z, Ta, Tb) = P

    # Negation in R2
    def R2neg(P):
        (N, D, E, F) = P
        return (D, N, E, GFp2.neg(F))

    # Constant-time point selection (in R2)
    def selectpt(c, P1, P2):
        (N1, D1, E1, F1) = P1
        (N2, D2, E2, F2) = P2
        return (
            GFp2.select(c, N1, N2),
            GFp2.select(c, D1, D2),
            GFp2.select(c, E1, E2),
            GFp2.select(c, F1, F2),
        )

    # Compute endomorphism images and negatives
    # XXX: Should check whether this table is for this point
    T = table
    if not T:
        T = table_endo(P)
    nT = [R2neg(P) for P in T]

    # Compute the scalar decomposition and recoding
    scalars = decompose(m)
    (s, d) = recode(scalars)

    # Compute the product
    Q = R2toR4(selectpt(s[64], T[d[64]], nT[d[64]]))
    for i in range(63, -1, -1):
        Q = DBL(Q)
        Q = ADD(Q, selectpt(s[i], T[d[i]], nT[d[i]]))

    return Q

########## Diffie-Hellman ##########

def DH_core(m, P, mul, table=None):
    if not PointOnCurve(P):
        raise Exception("Point not on curve")

    P0 = AffineToR1(P[0], P[1])
    P1 = DBL(P0)
    P2 = ADD(P1, R1toR2(P0))
    P3 = DBL(DBL(DBL(DBL(P2))))
    Q = ADD(P3, R1toR2(P0))
    Q = DBL(DBL(DBL(Q)))

    Q = R1toAffine(mul(m, Q, table=table))

    if Q == (Ox, Oy):
        raise Exception("DH computation resulted in neutral point")

    return Q

def DH_windowed(m, P, table=None):
    return DH_core(m, P, MUL_windowed, table=table)

def DH_endo(m, P, table=None):
    return DH_core(m, P, MUL_endo, table=table)

########## Self test for correctness ##########
########## (cases are from FourQlib) ##########

def test_definitions():
    test.test("0-on-curve", PointOnCurve((Ox, Oy)), True)
    test.test("G-on-curve", PointOnCurve((Gx, Gy)), True)

def test_encode():
    Genc = "87b2cb2b46a224b95a7820a19bee3f0e5c8b4c8444c3a74942020e63f84a1c6e"

    encTest = encode(Gx, Gy)
    test.test("encode", Genc, str(encTest).encode("hex"))

    try:
        decTest = decode(bytearray(Genc.decode("hex")))
        test.test("decode", (Gx, Gy), decTest)
    except Exception as err:
        test.test("decode", err, None)


def test_reps():
    x  = (0, 1)
    x2 = (0, 2)
    y  = (2, 0)
    y2 = (4, 0)
    xy = (2, 1)
    yx = (2, p1271 - 1)
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

    test.test("R1toR2", R1toR2(r1), r2)
    test.test("R1toR3", R1toR3(r1), r3)
    test.test("R2toR4", R2toR4(r2), r4)

def test_core():
    TEST_LOOPS = 1000

    # Test doubling
    A = (Gx, Gy, GFp2.one)
    for i in range(TEST_LOOPS):
        A = DBL(A)[:3]
    doubleP = ((0x2C3FD8822C82270FC9099C54855859D6, 0x4DA5B9E83AA7A1B2A7B3F6E2043E8E68),
               (0x2001EB3A576883963EE089F0EB49AA14, 0x0FFDB0D761421F501FEE5617A7E954CD))
    test.testpt("double", A, doubleP)

    # Test that the neutral element is neutral
    G = AffineToR1(Gx, Gy)
    O = AffineToR1(Ox, Oy)
    PP = ADD(G, R1toR2(O))
    test.testpt("neutral-r", PP, G)
    PP = ADD(O, R1toR2(G))
    test.testpt("neutral-l", PP, G)

    # Test point doubling by addition
    P = G
    for i in range(TEST_LOOPS):
        Q = R1toR2(P)
        P = ADD(P, Q)
    test.testpt("double-add", P, doubleP)

    # Test repeated addition of the same point
    P = G
    Q = R1toR2(P)
    P = DBL(P[:3])
    for i in range(TEST_LOOPS):
        P = ADD(P, Q)
    P1000 = ((0x3E243958590C4D906480B1EF0A151DB0, 0x5327AF7D84238CD0AA270F644A65D473),
             (0x3EF69A49CB7E02375E06003D73C43EB1, 0x293EB1E26DD23B4E4E752648AC2EF0AB))
    test.testpt("addition", P, P1000)

def test_mul(label, mul):
    TEST_LOOPS = 1000

    scalar = [0x3AD457AB55456230, 0x3A8B3C2C6FD86E0C, 0x7E38F7C9CFBB9166, 0x0028FD6CBDA458F0]
    coeff = range(TEST_LOOPS)
    for i in range(TEST_LOOPS):
         scalar[1] = scalar[2]
         scalar[2] += scalar[0]
         scalar[2] &= 0xffffffffffffffff
         coeff[i] = (scalar[0] << (0 * 64)) + (scalar[1] << (1 * 64)) + \
                    (scalar[2] << (2 * 64)) + (scalar[3] << (3 * 64))

    A = AffineToR1(Gx, Gy)
    for i in range(TEST_LOOPS):
        A = mul(coeff[i], A)

    mulP = ((0x257C122BBFC94A1BDFD2B477BD494BEF, 0x469BF80CB5B11F01769593547237C459),
            (0x0901B3817C0E936C281C5067996F3344, 0x570B948EACACE2104FE8C429915F1245))
    test.testpt(label, A, mulP)

def test_mul_windowed():
    # Test multiplication by one and two
    A = AffineToR1(Gx, Gy)
    B = MUL_windowed(1, A)
    test.testpt("mul-windowed-*1", B, A)
    A2 = DBL(A)
    B2 = MUL_windowed(2, A)
    test.testpt("mul-windowed-*2", B2, A2)

    # Test multiply over several iterations
    test_mul("mul-windowed", MUL_windowed)

    # Test fixed-based multiply
    T = table_windowed(A)
    B = MUL_windowed(1, A, table=T)
    B2 = MUL_windowed(2, A, table=T)
    test.testpt("mul-windowed-fixed-*1", B, A)
    test.testpt("mul-windowed-fixed-*2", B2, A2)

    failed = 0
    for i in range(10):
        m = getrandbits(256)
        B1 = MUL_windowed(m, A, table=T)
        B2 = MUL_windowed(m, A)
        if B1 != B2:
            failed += 1
    if failed == 0:
        print "[PASS] mul-windowed-fixed-rand"
    else:
        print "[FAIL] mul-windowed-fixed-rand"

def test_endo():
    TEST_LOOPS = 1000

    # Test phi endomorphism
    P = AffineToR1(Gx, Gy)
    for i in range(TEST_LOOPS):
        P = phi(P)
    phiP = ((0x5550AAB9E7A620EED5B5A3061287DB16, 0x3E61EBB9A1CB0210EC321E6CF33610FC),
            (0x5474BF8EC55603AE7E2851D5A8E83FB9, 0x5476093DBF8BF6BFA5077613491788D5))
    test.testpt("phi", P, phiP)

    # Test psi endomorphism
    P = AffineToR1(Gx, Gy)
    for i in range(TEST_LOOPS):
        P = psi(P)
    psiP = ((0x75AF54EDB41A2B93D8F3C8C24A2BC7E2, 0x065249F9EDE0C7984DE2466701F009A9),
            (0x06DBB85BFFB7C21E1C6E119ADD608104, 0x060A30903424BF13FD234D6C4CFA3EC1))
    test.testpt("psi", P, psiP)

def test_recoding():
    TEST_LOOPS = 1000

    # Test decomposition
    test_cases = [
        [0x92990788d66bf558052d112f5498111747b3e28c55984d43fed8c8822ad9f1a7L,
         [0xa8ea3f673f711e51L, 0xa08d1eae0b9e071dL, 0x55c8df690050276fL, 0x6396739dda88830f]],
        [0x48e5ca2a675ab49ca214b884813935024b0c61edc8d1305fe5230df341623348L,
         [0xa53ec4631945b875L, 0x521c0ba1261c1934L, 0x5c50ce912909185cL, 0x93b3c70960b44bad]],
        [0xae20e251c36cfa5be4d9f3d5a5edfed305a1e8f7f6394d9be58a15c4b0f1c5e9L,
         [0xa621ada9b3499c9fL, 0x7cd17e0095e7aae6L, 0x6e8d23b5bd10bb43L, 0x7f18c69f3025234c]],
        [0xb2c950abc87a55442cc00f1e3ac38f81b7e95036fd191ea134ff616d9806e10cL,
         [0x9b30a872ebea83afL, 0x8f6c73350447c9c3L, 0x72fdc76e3456d087L, 0x6ba39ba159b0c13d]],
        [0x8e2958a1475ed70762340e9797788e0061f21fcebd67889fdd4f4ce2b5f6b2deL,
         [0xbe8f3583a0934333L, 0xab45bf6d1bf80b37L, 0x4a19fc5cffe97809L, 0x5ea3baf1a1206442]],
    ]

    for t in test_cases:
        scalars = decompose(t[0])
        test.test("decompose", scalars, t[1])

    # Test recoding
    passed = 0
    failed = 0
    for n in range(TEST_LOOPS):
        k = getrandbits(256)
        scalars = decompose(k)
        (sign_masks, digits) = recode(scalars)

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
            failed += 1
        else:
            passed += 1
    if failed == 0:
        print "[PASS] recode"
    else:
        print "[FAIL] recode ({} / {})".format(failed, TEST_LOOPS)

def test_mul_endo():
    # Test multiplication by one and two
    A = AffineToR1(Gx, Gy)
    B = MUL_endo(1, A)
    test.testpt("mul-endo-*1", B, A)
    A2 = DBL(A)
    B2 = MUL_endo(2, A)
    test.testpt("mul-endo-*2", B2, A2)

    # Test multiply over several iterations
    test_mul("mul-endo", MUL_endo)

    # Test fixed-based multiply
    T = table_endo(A)
    B = MUL_endo(1, A, table=T)
    B2 = MUL_endo(2, A, table=T)
    test.testpt("mul-endo-fixed-*1", B, A)
    test.testpt("mul-endo-fixed-*2", B2, A2)

    failed = 0
    for i in range(10):
        m = getrandbits(256)
        B1 = MUL_endo(m, A, table=T)
        B2 = MUL_endo(m, A)
        if B1 != B2:
            failed += 1
    if failed == 0:
        print "[PASS] mul-windowed-fixed-rand"
    else:
        print "[FAIL] mul-windowed-fixed-rand"

def test_dh():
    TEST_LOOPS = 10

    def dhtest(label, dh):
        # Test that DH(m, P) == [392*m]P
        P = (Gx, Gy)
        failed = 0
        for i in range(TEST_LOOPS):
            m = getrandbits(256)
            Q1 = dh(m, P)
            Q2 = R1toAffine(MUL_windowed(392 * m, AffineToR1(P[0], P[1])))
            if Q1 != Q2:
                failed += 1
            P = Q1
        if failed == 0:
            print "[PASS] DH-{}-392".format(label)
        else:
            print "[FAIL] DH-{}-392".format(label)

        # Test that DH has the symmetry property
        G = (Gx, Gy)
        failed = 0
        for i in range(TEST_LOOPS):
            a = getrandbits(256)
            b = getrandbits(256)
            abG = dh(a, dh(b, G))
            baG = dh(b, dh(a, G))
            if abG != baG:
                failed += 1
        if failed == 0:
            print "[PASS] DH-{}-symm".format(label)
        else:
            print "[FAIL] DH-{}-symm {}".format(label, failed)

    dhtest("windowed", DH_windowed)
    dhtest("endo", DH_endo)

    def dhtest_fixed_base(label, dh, P, table):
        failed = 0
        A = AffineToR1(P[0], P[1])
        for i in range(TEST_LOOPS):
            m = getrandbits(256)
            Q1 = dh(m, P, table=table)
            Q2 = dh(m, P)
            if Q1 != Q2:
                failed += 1
        if failed == 0:
            print "[PASS] DH-{}-fixed".format(label)
        else:
            print "[FAIL] DH-{}-fixed {}".format(label, failed)

    G = (Gx, Gy)
    G392 = MUL_endo(392, AffineToR1(Gx, Gy))
    T_windowed = table_windowed(G392)
    T_endo = table_endo(G392)
    dhtest_fixed_base("windowed", DH_windowed, G, T_windowed)
    dhtest_fixed_base("endo", DH_endo, G, T_endo)

    # Test that points not on the curve are rejected
    try:
        DH_endo(1, ((0, 0), (0, 0)))
        print "[FAIL] DH-reject-not-on-curve"
    except:
        print "[PASS] DH-reject-not-on-curve"

    # Test that 392-torsion points are rejected
    P392 = ((0x1318020702de23bc3c9b73c751b4b192, 0x77ab39a7d8990c0a18e3c409fbd81a95),
            (0x515854b6d19cc2da1ea2b43b5121a22e, 0x763f89e129497361d74dff5063e66682))
    try:
        DH_endo(1, P392)
        print "[FAIL] DH-reject-392-torsion"
    except:
        print "[PASS] DH-reject-392-torsion"

if __name__ == "__main__":
    import test
    test_definitions()
    test_encode()
    test_reps()
    test_core()
    test_mul_windowed()
    test_endo()
    test_recoding()
    test_mul_endo()
    test_dh()
