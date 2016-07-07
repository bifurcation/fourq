---
title: Curve4Q
abbrev: 
docname: draft-barnes-cfrg-4q
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
 -
       ins: W. Ladd
       name: Watson Ladd
       organization: UC Berkeley
       email: watsonbladd@gmail.com

informative:
   EFD:
      target: https://hyperelliptic.org/EFD/
      title: Explicit-Formulas Database
      author:
        -
                ins: D.J. Bernstein
        -
                ins: T. Lange
   Curve4Q:
      target: https://eprint.iacr.org/2015/565.pdf
      title: "FourQ:four-dimensional decompositions on a Q-curve over the Mersenne prime"
      author:
        -
                ins: C. Costello
        -
                ins: P. Longa
                
   TwistedRevisted:
      target: http://iacr.org/archive/asiacrypt2008/53500329/53500329.pdf
      title: Twisted Edwards Curves Revisted
      author:
        -
            ins: H. Hisil
        -
            ins: K.H.Wong
        -
            ins: G. Carter
        -
            ins: E. Dawson
      
--- abstract

This document specifies an elliptic curve over a quadratic extension
of a prime field that offers the fastest known Diffie-Hellman key
agreements while using compact keys. This high performance does not
require vectorization and applies to signature verification. The best
known attacks require 2^125 or so operations, comparable to X25519, while
performance is twice as good.

--- middle

# Introduction

Public key cryptography continues to be computationally expensive particularly on less powerful devices. While
recent advances have substantially reduced the cost of elliptic curve operations, the use of endomorphisms enables
even more speedups.

## Comparison to prior X25519 / X448

* Similar security level, also constant-time
* Faster due to endomorphisms
* More magic: potentially risky.

# Four-dimensional decompsitions

As described in [Curve4Q], the elliptic curve described in this document is
the only known curve that permits a four dimensional GLS-GLV decomposition
over the highly efficient field GF(p^2) with p=2^127-1.

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

Addition of two elements A=a0+a1*i, B=b0+b1*i is performed coordinatewise: A+B=(a0+b0)+(a1+b1)*i.
Multiplication is similarly simple: A*B=(a0b0-a1b1)+(a0b1+a1b0)*i. Lastly there is a field automorphism
conj(A)=a0-a1*i.

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

This algorithm takes a scalar m and a point P which is a N torsion point and computes m*P.
It computes phi(P), psi(P), and
psi(phi(P)), where phi and psi are endomorphisms and takes
m*P=a_0*P+a_1*phi(P)+a_2*psi(P)+a_3*psi(phi(P)), where a_0, a_1,
a_2, and a_3 have been computed as below. It is far more efficient then the above one.
In its description we make use of constants listed in an appendix.

### Alternative Point Representations and addition laws

We use the following 3 representations of a point (x, y) on the
curve. All representations use X, Y, Z satisfying x=X/Z, y=Y/Z. The
point at infinity is (1,1,0). These representations differ in auxiliary
data used to speed some operations. By omitting their computation when
they are not needed we save operations. In three of these representations
T=XY/Z is used to define them.

R1 points are (X,Y,Z,Ta,Tb) where T=Ta*Tb. R2 points are
(N, D, E, F)=(X+Y,Y-Z,2Z,2dT). R3 are (N, D, Z, T)=(X+Y,Y-X,Z,T), and R4 is (X,Y,Z).  A point
doubling takes an R4 point, and produces an R1 point.  There are two
kinds of addition: ADD_core eats an R2 and an R3 point and produces an
R1 point, and ADD eats an R1 and an R2 point, by converting the R1
point into an R3 point and then proceeding. Exposing these two
operations and the multiple representations helps save time in
precomputing tables and in using the tables effectively. The conversions
between point formats are obvious.

These operations have the following explicit formulas developed by
many people over the years ([EFD],[TwistedRevisted]). We present
the operations as functions in pseudocode.

Doubling is computed as follows

~~~~
DBL(X1, Y1, Z1):
  A = X1^2
  B = Y1^2
  C = 2*Z1^2
  D = -1*A
  E = (X1+Y1)^2-A-B
  G = D+B
  F = G-C
  H = D-B
  X3 = E*F
  Y3 = G*H
  Ta3 = E
  Tb3 = H
  Z3 = F*G
return(X3, Y3, Z3, Ta3, Tb3)
~~~~

ADD_core is computed as follows:

~~~~
ADD_core(N1, D1, E1, F1, N2, D2, Z2, T2):
   A = D1*D2
   B = N1*N2
   C = T2*F1
   D = Z2*E1
   E = B-A
   F = D-C
   G = D+C
   H = B+A
   X3 = E*F
   Y3 = G*H
   Ta3 = E
   Tb3 = H
   Z3 = F*G
return (X3, Y3, Z3, Ta3, Tb3)
~~~~

### Endomorphisms and Isogenies

Our endomorphisms and isogenies mostly work in projective coordinates. We present formulas for
tau and tau_dual, and then for upsilon and chi. phi(Q)=tau_dual(upsilon(tau(Q)), where of course
tau is performed first and similarly psi(Q)=tau_dual(chi(tau(Q))). tau_dual outputs points in R1,
tau consumes affine points, and chi and upsilon input and output points in projective coordinates on a different curve.

~~~~
tau(X1, Y1):
   A = X1^2
   B = Y1^2
   C = A+B
   D = A-B
   X2 = ctau1*X1*Y1*D
   Y2 = -(2*Z1^2+D)*C
   Z2 = C*D
return(X2, Y2, Z2)
~~~~

~~~~
 tau_dual(X1, Y1, Z1):
  A = X1^2
  B = Y1^2
  C = A+B
  T2a = B-A
  D = 2*Z1^2-Ta
  T2b = ctaudual1*X1*Y1
  X2 = T2b*C
  Y2 = D*T2a
  Z2 = D*C
return(X2, Y2, Z2)
~~~~

~~~~
upsilon(X1, Y1, Z1):
   A = cphi0*X1*Y1
   B = Y1*Z1
   C = Y1^2
   D = Z1^2
   F = D^2
   G = B^2
   H = C^2
   I = cphi1*B
   J = C+cphi2*D
   K = cphi8*G+H+cphi9*F
   X2 = conj(A*K*(I+J)*(I-J))
   L = C+cphi4*D
   M = cphi3*B
   N = (L+M)*(L-M)
   Y2 = conj(cphi5*D*N*(H+cphi6*G+cphi7*F))
   Z2 = conj(B*K*N)                                      
return(X2, Y2, Z2)
~~~~

~~~~
chi(X1, Y1, Z1):
   A = conj(X1)  
   B = conj(Y1) 
   C = conj(Z1)^2   
   D = A^2 
   F = B^2 
   G = B*(D+cpsi2*C)
   H = -(D+cpsi4*C)
   X2 = cpsi1*A*C*H
   Y2 = G*(D+cpsi3*C)
   Z2 = G*H 
return(X2, Y2, Z2)
~~~~



### Scalar Recoding

Scalar recoding has two parts. The first is to decompose the scalar into four small integers, the second
is to encode these integers into an equivalent form satisfying certain properties, which will be what is used
by the multiplication algorithm.

This decomposition uses another bunch of constants defining four
vectors with integer coordinates b1, b2, b3, b4. These constants are
64 bit. Then there are four constants l1, l2, l3, l4 which are long
integers used to implement rounding.


Let c=2b1-b2+5b3+2b4 and c'=2b1-b2+5b3+b4. Then compute ti=floor(li*m/2^256), and then compute
a=(a1, a2, a3, a4)=(m,0,0,0)-t1*b1-t2*b2-t3*b3-t4*b4. Precisely one of a+c and a+c' has an odd first
coordinate: this is the one fed into the next scalar recoding step. Each entry is 64 bits after this
calculation, and so the ti and m can be truncated to 64 bits.

The second step takes the four 64 bit integers a1, a2, a3, a4
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

Next comes table initialization. Let Q=psi(P), R=phi(P), S=psi(phi(P)), all in R2, and take P in
R3 form.

T will be a table of 8 R2 format points.
T[0] is P in R2
T[1] is P+Q
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

We note that taking the inverse of a point is simply negating the y coordinate.
This multiplication algorithm is only valid for N torsion points. It will produce the
wrong answer for points that are not N torsion points, and in the process cause massive
damage to security. Implementations MUST NOT use this algorithm on anything that is not
a torsion point.

# Use of the scalar multiplication primitive in Diffie-Hellman

The above scalar multiplication primitive can be used to implement elliptic curve Diffie-Hellman
with cofactor. To multiply by the cofactor of 392 requires nine doublings and two additions.

~~~
DH(m, P):
      Check P is on the curve: if not return failure
      Q = 392*P
      Compute m*Q with the multiplication algorithm
return m*Q
~~~~

Two users, Alice and Bob, can carry out the following steps to derive a shared key:
Both pick a random string of 32 bytes, mA and mB respectively. Alice computes
A = DH(mA, G), Bob B = DH(mB, G). They exchange A and B, and then Alice computes
KAB = DH(mA, B), while Bob computes KBA = DH(mB, A). The coordinates of G are
found in the appendix. [[TODO: test vector and make this true]]

It is possible to have an even more efficient fixed-base multiplication, either by storing
the table that the above routine uses in memory, or via comb-based methods.

# IANA Considerations

IANA need take no action.

# Security Considerations

Claus Diem has steadily reduced the security of elliptic curves
defined over extension fields of degree greater then two over large
characteristic fields. There is considerable concern about curves like
FourQ defined over extension fields, which are considered less
conservative then curves over prime fields. While the best attack for
Diffie-Hellman on this curve remains generic, this may change.

It is absolutely essential that points input to scalar multiplication
algorithms are checked for being on the curve first. Removing such
checks may result in revealing the entire scalar to an attacker. Implementations
MUST check that the points are on the curve. The curve is not twist-secure: single
coordinate ladders MUST NOT operate on the twist but MUST validate points before
operating on them.

The arithmetic operations and table loads must be done in constant
time to prevent timing attacks. Side-channel analysis is a constantly
moving field, and implementors must be extremely careful.

--- back

# Constants
ctau1= 221360928884514619410*i + 33754435779700894835198039471158097091

ctaudual1 = 170141183460469231510326374831369486353*i + 99231301967130569661901792840482943028

cphi0 = 49615650983565284830950896420241471514*i + 110680464442257309687

cphi1 = 131306912742858181648727312260439119609*i + 92233720368547758087

cphi2 = 160666015865631300014011952927357137809*i + 276701161105643274261

cphi3 = 107027644557995218531204623577807990436*i + 36893488147419103235

cphi4 = 24279268184862963117522688682631129173*i + 55340232221128654851

cphi5 = 92472642025247131565767320804994133491*i + 184467440737095516175

cphi6 = 14804100590025031399847337894104161255*i + 332041393326771929112

cphi7 = 76283848507754718862858058709728786458*i + 442721857769029238819

cphi8 = 41635071732389019719735756359456329456*i + 3135946492530623774960

cphi9 = 21045324596686230484035983431638590725*i + 39844967199212631493615

cpsi1 = 4095177184363520459066*i + 57123674603396429897431647433607300847

cpsi2 = 44824135016688633386011024159913800562*i + 4205857648805777768771

cpsi3 = 101947809620085063283442671593521101409*i + 110680464442257309705

cpsi4 = 68193373840384168448244632122363004318*i + 170141183460469231621006839273626796022

l1 = 50127518246259276682880317011538934615153226543083896339791

l2 = 22358026531042503310338016640572204942053343837521088510715

l3 = 5105580562119000402467322500999592531749084507000101675068

l4 = 19494034873545274265741574254707851381713530791194721254848

b1 = [650487742939046294, -1397215820276968864, 523086274270593807,
   -598824378691085905]

b2 = [2110318963211420372, -1, 1, 2727991412926801872]

b3 = [1705647224544756482, 199320682881407569,
   -3336360048424633503, 765171327772315031]

b4 = [1400113754146392127, 3540637644719456050, -471270406870313397, -1789345740969872106]
