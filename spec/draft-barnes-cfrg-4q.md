---
title: Curve4Q
abbrev: 
docname: draft-rescorla-tls-subcerts-latest
category: std

ipr: trust200902
area: Security
workgroup: 
keyword: Internet-Draft

stand_alone: yes
pi: [toc, sortrefs, symrefs]

author:
 -
       ins: R. Barnes
       name: Richard Barnes
       organization: Mozilla
       email: rlb@ipv.sx


--- abstract

This document specifies an elliptic curve over a quadratic extension
of a prime field that offers the fastest known Diffie-Hellman key
agreements while using compact keys. This high performance does not
require vectorization and applies to signature verification. The best
known attacks require 2^125 or so operations, comparable to X25519.

--- middle

# Introduction

[[ Ed. - TODO ]]

## Comparison to prior X25519 / X448

* Similar security level, also constant-time
* Faster due to endomorphisms
* More magic: potentially risky.

# Four-dimensional decompsitions

As described in [Curve4Q], the elliptic curve described in this document is
the only known curve that satisfies the following requirements:

* GF(p^2) for p = 2^127-1
* 4-dimensional decomposition
* Security of around 2^120 operations for discrete log computation

# Mathematical Prerequisites

Curve4Q is defined over the finite field GF(p^2), where p is the Mersenne prime
2^127 - 1.  Elements of this finite field have the form (a + b * i), where a and
b are elements of the finite field GF(p) (i.e., integers mod p) and i^2 = -1.

Curve4Q is the twisted Edwards curve over GF(p^2) defined by the following curve
equation:

~~~~~
Y^2 - X^2 = 1 + d * X^2 * Y^2, with

d = 0x00000000000000e40000000000000142 +
    0x5e472f846657e0fcb3821488f1fc0c8d * i
~~~~~

This order of this curve is 23 ·72 ·N, where N is the following 246-bit prime:

~~~~~
N = 0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe0f79992fb2540ec7768ce7
~~~~~


# Curve Points

Elements of GF(p) are represented as 16-octet little-endian integers.  An element
x0 + x1*i of GF(p^2) is represented on the wire by the concatenation of the byte strings for
x0 and x1. Implementations will use whatever internal representation they desire, but we will
describe the operations on elements of GF(p^2) assuming that X=x0+x1*i, with each x0 and x1
an element of GF(p).

A point on this curve is serialized as a sequence of 64 octets, representing
the point's coordinates X = x0 + x1*i and Y = y0 + y1*i.  The individual
coordinates are serialized as little-endian integers.

~~~~~
|--------------- X ---------------|--------------- Y ---------------|
|        x0      |        x1      |        y0      |        y1      |
|................|................|................|................|
~~~~~

[[ Ed. - Endianness? ]]

[[ Ed. - Can we get by with just X? ]]
[[WBL: Nope. p^2-1 is divisible by 2^128. This makes taking a square root Lots Of Fun.]]

Addition of two elements A=a0+a1*i, B=b0+b1*i is performed coordinatewise: A+B=(a0+b0)+(a1+b1)*i.
Multiplication is similarly simple: A*B=(a0b0-a1b1)+(a0b1+a1b0)*i. Lastly there is a field automorphism
Frob_p(A)=a0-a1*i.

In the Python code samples below, we represent elements of GF(p^2) as Python
tuples, with two elements, (x0, x1) = x0 + x1*i.  Likewise, points are
represented by tuples of field elements (X, Y).

~~~~~
<CODE BEGINS>
def decodeLittleEndian(b, bits):
    return sum([b[i] << 8*i for i in range((bits+7)/8)])

def encodeLittleEndian(b, bits):
    return bytearray([(x >> (8*i)) & 0xff for i in range(bits >> 3)])

def decodeGFp2(b):
    x0 = decodeLittleEndian(b, 128)
    x1 = decodeLittleEndian(b[16:], 128)
    return (x0, x1)

def encodeGFp2(X):
    b0 = encodeLittleEndian(X[0], 128)
    b1 = encodeLittleEndian(X[1], 128)
    return b0 + b1

def decodePoint(b):
    X = decodeGFp2(b)
    Y = decodeGFp2(b[32:])
    return (X, Y)

def encodePoint(P):
    B0 = encodeGFp2(P[0])
    B1 = encodeGFp2(P[1])
    return B0 + B1
<CODE ENDS>
~~~~~


# The Curve4Q Function

The Curve4Q function produces a 64-octet string from an input point P and a
256-bit integer coefficient m.  The output of the Curve4Q function
are the coordinates of the curve point m*P, encoded as described above.

~~~~~
Curve4Q(m, P) = encodeGFp2(MUL(m, P)[0])
~~~~~

The function encodeGFp2 is defined above.  The MUL function represents scalar
multiplication according to the group law of the curve.  We give two explicit
algorithms for computing MUL below: A baseline algorithm that is short, simple,
and slow, and an optimized algorithm that uses some precomputation to achieve
significant speed benefits.


## Baseline Point Multiplication Algorithm

The function defined here represents a constant-time implementation of
"textbook" scalar multiplication on the curve.  It is presented mainly as a
reference for implementers who might want a simpler implementation. The algorithm
in the next section provides substantially greater performance.

This code uses formulas from [TwistedRevisted].
~~~~~
Inputs:
- A curve point P = (X, Y)
- A 256-bit integer m

Pxy = X + Y
Pyx = Y - X
P2z = 2
P2dt = 2 * d * X * Y

Sx = 0
Sy = 1
Sz = 1

for t = 255 down to 0:
  m_t = (m >> t) & 1
  // Constant-time selection; see below
  Axy = cselect(m_t, Pxy, 1)
  Ayx = cselect(m_t, Pyx, 1)
  A2z = cselect(m_t, P2z, 2)
  A2dt = cselect(m_t, P2dt, 0)

  XX = Sx * Sx
  YY = Sy * Sy
  ZZ = Sz * Sz

  A = Sx + Sy
  AA = A * A
  B = XX + YY
  C = YY - XX
  D = AA - B
  E = 2*ZZ - C
  F = D * E
  G = B * C

  Txy = F + G
  Tyx = G - F
  Tz = C * E
  Tt = D * B

  TAxy = Txy * Axy
  TAyx = Tyx * Ayx
  TA2z = T2z * A2z
  TAt = Tt * A2dt

  H = TAz + TAt
  I = TAz - TAt
  K = TAxy + TAyx
  L = TAxy - TAyx

  Sx = L * I
  Sy = K * H
  Sz = H * I

Sz = 1 / Sz
Sx = Sx * Sz
Sy = Sy * Sz

Return (Sx, Sy)
~~~~~

The cselect function returns its second or third argument depending on whether
the first argument is one or zero, respectively.  This function SHOULD be
implemented in constant time.  For example, this can be done as follows:

~~~~~
cselect(c, x, y) = ((x ^ y) & mask(c)) ^ y
~~~~~

Where mask(c) is the all-1 or all-0 word of the same length as x and y,
computed, e.g., as mask(c) = 0 - c.

## Optimized Point Multiplication Algorithm

This algorithm takes a scalar m and a point P and computes 392*m*P. As
392 is the cofactor, it suffices to check that P is on the curve to
prevent small-subgroup attacks. This algorithm uses isogenies to speed
up the scalar multiplication. The algorithm computes Q=392*P via a
simple addition chain. It then computes \phi(Q), \psi(Q), and
\psi(\phi(Q)), where \phi and \psi are endomorphisms and takes
mQ=a_0*P+a_1*\phi(Q)+a_2*\psi(Q)+a_3*\psi(\phi(Q)), where a_0, a_1,
a_2, and a_3 have been computed as below. It is far more efficient then the above one

### Alternative Point Representations and addition laws

We use the following 3 representations of a point (x, y) on the
curve. All representations use X, Y, Z satisfying x=X/Z, y=Y/Z. The
point at infinity is (1,1,0). These representations differ in auxiliary
data used to speed some operations. By omitting their computation when
they are not needed we save operations. In three of these representations
T=XY/Z is used to define them.

R1 points are (X,Y,Z,Ta,Tb) where T=Ta*Tb. R2 points are
(X+Y,Y-Z,2Z,2dT). R3 are (X+Y,Y-X,Z,T), and R4 is (X,Y,Z).  A point
doubling takes an R4 point, and produces an R1 point.  There are two
kinds of addition: ADD_core eats an R2 and an R3 point and produces an
R1 point, and ADD eats an R1 and an R2 point, by converting the R1
point into an R3 point and then proceeding. Exposing these two
operations and the multiple representations helps save time in
precomputing tables and in using the tables effectively.

These operations have the following explicit formulas largely taken from [EFD]
[TODO]

### Endomorphisms and Isogenies
Our endomorphisms and isogenies mostly work in projective coordinates. We present formulas for
\tau and \hat{\tau}, and then for \upsilon and \chi. \phi=\hat{\tau}\upsilon\tau, where of course
\tau is performed first and similarly \psi=\hat{\tau}\chi\tau. \hat{\tau} outputs points in R1 while
\chi and \upsilon input and output points in projective form

We begin by defining some constants [WBL: figure out how represented and convert to human form]

Each of our formulas consumes X1, Y1, Z1 and produces X2, Y2, Z2  or X2, Y2, Z2, T2a, T2b. Note that
the outputs are not on the curve: only \psi and \phi produce outputs on the curve.

The following is the formula for \tau:
    X2 = ctau*X1*Y1*(X1^2-Y1^2)
    Y2 = (X1^2+Y1^2)*(-2*Z1^2-(X1^2-Y1^2))
    Z2 = (X1^2+Y1^2)*(X1^2-Y1^2)

\tau is most efficiently computed by
    A = X1*Y1
    B = (X1+Y1)^2-2A
    C = (X1-Y1)^2-2A
    X2 = A*C
    Y2 = A*(-2Z1^2-C)
    Z2 = B*C

[TODO] More here
### Scalar Recoding

Scalar recoding has two parts. The first is to decompose the scalar into four small integers, the second
is to encode these integers into an equivalent form satisfying certain properties, which will be what is used
by the multiplication algorithm.

This decomposition uses another bunch of constants defining four
vectors with integer coordinates b1, b2, b3, b4. These constants are
64 bit. Then there are four constants l1, l2, l3, l4 which are long
integers used to implement rounding.

[TODO: determine why the ti can be 64 bits. Possible magic]
Let c=2b1-b2+5b3+2b4 and c'=2b1-b2+5b3+b4. Then compute ti=floor(li*m/2^256), and then compute
a=(a1, a2, a3, a4)=(m,0,0,0)-t1*b1-t2*b2-t3*b3-t4*b4. Precisely one of a+c and a+c' has an odd first
coordinate: this is the one fed into the next scalar recoding step.

The scalar recoding step takes the four 64 bit integers a1, a2, a3, a4
from the previous step and outputs two arrays m[0]..m[64] and
d[0]..d[64]. Each entry of d is between 0 and 7, and each entry in m
is -1 or 0. Rather then describe the properties required of this
encoding, we present the algorithm. bit(x, n) denotes the nth bit of x.
~~~~~
m[64]=-1
for i=0 to 63 do
    d[i] = 0
    m[i] = -bit(a1, i+1)
    for j = 2 to 4 do:
    	d[i] = d[i]+bit(aj, 0)<<(j-2)
	c = (bit(a1, i+1)|bit(aj,0)) xor bit(a1, i+1)
	aj = aj/2+c
d[64]=a2+2a3+4a
~~~~~
### Multiplication
We begin by taking the input point P and multiplying by 392 using a double and add routine.

Next we recode the scalar m into d[i] and m[i].

Next comes table initialization. Let Q=\psi(P), R=\phi(P), S=\psi(\phi(P)), all in R2, and take P in
R3 form.

T will be a table of 8 R2 format points.
T[0] is P in R2
T[1] is P+Q, again in R2
T[2] is R+P
T[3] is R+P+Q
T[4] is S+P
T[5] is S+P+Q
T[6] is S+P+R
T[7] is S+P+Q+R

By converting the table entries as soon as they are computed ADD_core can be used without need for
extraneous arithmetic operations.

Define s[i] to be 1 if m[i] is -1 and -1 if m[i] is 0. Then our multiplication algorithm is the following:
~~~~~
Q = s[64]*T[d[64]] in R4
for i=63 to 0 do:
    Q=DBL(Q)
    Q=ADD(Q, s[i]*T[di])
return Q
~~~~~
This multiplication algorithm has a regular pattern of operations and no exceptional cases.
[ TODO ]

# IANA Considerations

# Security Considerations

Claus Diem has steadily reduced the security
of elliptic curves defined over extension fields of odd prime fields
along with Sameav. There is considerable concern about curves like
FourQ defined over extension fields, which are considered less
conservative then curves over prime fields. While the best attack for
Diffie-Hellman on this curve remains generic, this may change.

It is absolutely
essential that points input to scalar multiplication algorithms are
checked for being on the curve first. Removing such checks may result
in revealing the entire scalar to an attacker.

The arithmetic operations and table loads must be done in constant
time to prevent timing attacks. Side-channel analysis is a constantly
moving field.

--- back
