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

This document specifies an elliptic curve over a quadratic extension of a prime
field that offers comparable security guarantees to existing curves, while
offering faster curve operations.  This curve is intended to operate at the
~128-bit security level, and is the unique curve satisfying a list of required
properties.

--- middle

# Introduction

[[ Ed. - TODO ]]

## Comparison to prior X25519 / X448

* Similar security level, also constant-time
* Faster due to endomorphisms
* Less magic

# Uniqueness

As described in [Curve4Q], the elliptic curve described in this document is
the only possible elliptic curve that satisfies the following requirements:

* GF(p^2) for p = 2^127-1
* 4-dimensional decomposition

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
x0 + x1*i of GF(p^2) is represented by the concatenation of the byte strings for
x0 and x1.

A point on this curve is represented as a sequence of 64 octets, representing
the point's coordinates X = x0 + x1*i and Y = y0 + y1*i.  The individual
coordinates are represented as little-endian integers.

~~~~~
|--------------- X ---------------|--------------- Y ---------------|
|        x0      |        x1      |        y0      |        y1      |
|................|................|................|................|
~~~~~

[[ Ed. - Endianness? ]]

[[ Ed. - Can we get by with just X? ]]

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
256-bit integer coefficient m.  The output of the Curve4Q function is the
X-coordinate of the curve point m*P, encoded as described above.

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
reference for implementers who might want a simpler implementation.  The
optimized version presented in the next section is much faster, and is the
recommended implementation.

[[ TODO: Something like Section 5 of RFC 7748 ]]

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

### Alternative Point Representations

[[ TODO ]]

### Scalar Recoding

[[ TODO ]]

### Endomorphisms

[[ TODO ]]

### Multiplication

[[ TODO ]]

# IANA Considerations

# Security Considerations

--- back
