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
 -
       ins: P. Longa
       name: Patrick Longa
       organization: Microsoft Research
       email: plonga@microsoft.com

informative:

    Curve4Q:
      target: "https://eprint.iacr.org/2015/565.pdf"
      title: "FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime"
      date: 2016
      author:
         -
              ins: C. Costello
         -
              ins: P. Longa

    Distinguished:
       target: "http://people.scs.carleton.ca/~paulv/papers/JoC97.pdf"
       title: "Parallel Collision Search with Cryptanalytic Applications"
       date: 1996
       author:
          -
             ins: P.C. van Oorschot
          -
             ins: M.J. Wiener

    Exceptional:
       target: "https://www.iacr.org/archive/pkc2003/25670224/25670224.pdf"
       title: "Exceptional procedure attack on elliptic curve cryptosystems"
       date: 2003
       author:
          -
              ins: T. Izu
          -
              ins: T. Takagi

    Invsqr:
        target: "http://eprint.iacr.org/2012/309.pdf"
        title: Fast and compact elliptic-curve cryptography
        date: 2012
        author:
           -
              ins: M. Hamburg

    FourQlib:
      target: "https://www.microsoft.com/en-us/research/project/fourqlib/"
      title: "FourQlib"
      date: 2016
      author:
         -
              ins: C. Costello
         -
              ins: P. Longa

    GLV:
       target: "https://www.iacr.org/archive/crypto2001/21390189.pdf"
       title: "Faster Point Multiplication on Elliptic Curves with Efficient Endomorphisms"
       date: 2001
       author:
          -
             ins: R. P. Gallant
          -
             ins: R. J. Lambert
          -
             ins: S. A. Vanstone

    GLS:
        target: "https://www.iacr.org/archive/eurocrypt2009/54790519/54790519.pdf"
        title: "Endomorphisms for Faster Elliptic Curve Cryptography on a Large Class of Curves"
        date: 2009
        author:
           -
             ins: S. D. Galbraith
           -
             ins: X. Lin
           -
             ins: M. Scott

    SQRT:
        target: "https://eprint.iacr.org/2012/685.pdf"
        title: "Square Root Computation over Even Extension Fields"
        date: 2012
        author:
           -
              ins: G. Adj
           -
              ins: F. Rodriguez-Henriquez

    SchnorrQ:
       target: "https://www.microsoft.com/en-us/research/wp-content/uploads/2016/07/SchnorrQ.pdf"
       title: "SchnorrQ"
       date: 2016
       author:
          -
             ins: C. Costello
          -
             ins: P. Longa

    TwistedRevisited:
       target: "http://iacr.org/archive/asiacrypt2008/53500329/53500329.pdf"
       title: Twisted Edwards Curves Revisited
       date: 2008
       author:
          -
              ins: H. Hisil
          -
              ins: K-H. Wong
          -
              ins: G. Carter
          -
              ins: E. Dawson

    Twisted:
       target: "http://eprint.iacr.org/2008/013.pdf"
       title: Twisted Edwards Curves
       date: 2008
       author:
          -
              ins: D. J. Bernstein
          -
              ins: P. Birkner
          -
              ins: M. Joye
          -
              ins: T. Lange
          -
              ins: C. Peters

--- abstract

This document specifies a twisted Edwards curve that takes advantage
of arithmetic over the field GF(2^127-1) and two endomorphisms to achieve
the speediest Diffie-Hellman key agreements over a group of order
approximately 2^246, which provides around 128 bits of security.
Curve4Q implementations are roughly twice as fast as those of Curve25519,
and when not using endomorphisms take eighty percent of the time taken by Curve25519.

--- middle

# Introduction

Public key cryptography continues to be computationally expensive, particularly
on less powerful devices. While recent advances in efficient formulas for
addition and doubling have substantially reduced the cost of elliptic curve
operations in terms of field operations, the number of group operations involved
in scalar multiplication has not been reduced in the curves considered for IETF
use. Using curves with an efficiently computable endomorphism reduces the number
of group operations by turning one long scalar multiplication into the sum of
two multiplications by scalars of half the length, which can be evaluated more
efficiently.

By using curves over quadratic extensions there are more endomorphism families
to pick from, and the field operations become more efficient compared to prime
fields of the same size.  The ideal case is given by curves equipped with two
distinct endomorphisms, so that it becomes possible to divide scalars into four
parts.  We also focus on curves defined over the field GF(p^2) for the Mersenne
prime p = 2^127 - 1, which offers extremely efficient arithmetic.  Together,
these improvements substantially reduce computation time compared to other
proposed Diffie-Hellman key exchange and digital signature schemes.  However,
the combined availability of these features severely restricts the curves that
can be used for cryptographic applications.

As described in {{Curve4Q}}, Curve4Q is the only known elliptic curve that
permits a four dimensional decomposition over the highly efficient field GF(p^2)
with p = 2^127 - 1 and has a prime order subgroup of order approximately 2^246,
which provides around 128 bits of security. No other known elliptic curve with
such a decomposition has a larger prime order subgroup over this field.  This
"uniqueness" allays concerns about selecting curves vulnerable to undisclosed
attacks.

Curve4Q can be used to implement Diffie-Hellman key exchange, as described
below.  It is also possible to use Curve4Q as the basis for digital signature
scheme (e.g., {{SchnoorQ}}).

# Mathematical Prerequisites

Curve4Q is defined over the finite field GF(p^2), where p is the Mersenne prime
2^127 - 1.  Elements of this finite field have the form (a + b * i), where a and
b are elements of the finite field GF(p) (i.e., integers mod p) and i^2 = -1.

Let A = a0 + a1\*i and B = b0 + b1\*i be two elements of GF(p^2). Below we
present formulas for computing addition, subtraction, multiplication, squaring,
conjugation and inversion.

~~~~
A + B = (a0 + b0) + (a1 + b1)*i

A - B = (a0 - b0) + (a1 - b1)*i

A * B = (a0*b0 - a1*b1) + ((a0+a1)*(b0+b1)-(a0*b0 - a1*b1))*i
      = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)*i

A * A = (a0+a1)*(a0-a1) + 2*a0*a1*i

conj(A) = a0 - a1*i

1/A = conj(A)/(a0^2+a1^2)
~~~~

Curve4Q is the twisted Edwards curve E over GF(p^2) defined by the
following curve equation:

~~~~~
E: -x^2 + y^2 = 1 + d * x^2 * y^2, with

d = 0x00000000000000e40000000000000142 +
0x5e472f846657e0fcb3821488f1fc0c8d * i
~~~~~

Let E(GF(p^2)) be the set of pairs (x, y) of elements of GF(p^2)
satisfying this equation. This set forms a group with the addition operation
shown below:
~~~~
                    /   x1*y2+y1*x2         y1*y2+x1*x2   \
(x1, y1)+(x2, y2) = |------------------,------------------|
                    \ 1+d*x1*y1*x2*y2    1-d*x1*y1*x2*y2  /

~~~~

As d is not a square in GF(p^2), and -1 is, this formula never involves a
division by zero when applied to points on the curve. That is, the formula is
*complete* and works without exceptions for any input in E(GF(p^2)).  The
identity element is (0, 1), and the inverse of (x,y) is (-x, y).  The order of
this group is \#E = 2^3 · 7^2 · N, where N is the following 246-bit prime:

~~~~~
N = 0x29cbc14e5e0a72f05397829cbc14e5dfbd004dfe0f79992fb2540ec7768ce7
~~~~~

Points P on E such that N\*P = (0, 1) are N-torsion points. Given a point P and
Q which are both N-torsion points, it is difficult to find m such that Q = m\*P.
This is the elliptic curve discrete logarithm problem, which is closely related
to the security of Diffie-Hellman key exchanges as the best known attacks on the
Diffie-Hellman problem involve solving the discrete logarithm problem. The best
known algorithms take approximately 2^123 group operations.

This group has two different efficiently computable endomorphisms, as described
in {{Curve4Q}}. As discussed in {{GLV}} and {{GLS}} these endomorphisms
multiplication by a large scalar to be computed using multiple multiplications
by smaller scalars, which can be evaluated in much less time overall.

# Representation of Curve Points

Elements a in GF(p) are represented as 16 byte little endian integers which are
the numbers in the range [0, p). The 16 bytes b[0], b[1],... b[15] represent
b[0]+256\*b[1]+256^2\*b[2]+...+256^15\*b[15].  Since we are representing numbers
in the range [0, 2^127-1), the top bit of b[15] is always zero.

An element x0 + x1\*i of GF(p^2) is represented on the wire by the concatenation
of the encodings for x0 and x1. A point (x, y) on Curve4Q is serialized in a
compressed form as the representation of y with a modified top bit. This top bit
is used to disambiguate between x and -x.

To carry out this disambiguation we order the elements of GF(p^2) as follows: to
compare x = x0+x1\*i with y = y0+y1\*i assuming all coordinates are in [0, p) we
compare x0 with y0, and, if they are equal, compare x1 with y1. This is the
lexicographic ordering on (x0, x1) where each value is in the range [0, p).

The high bit of a compressed point is 0 if the smaller possible x value is
correct, and 1 if the larger possible x value is correct.

~~~~~
|--------------- y ---------------|
|       y0     |0|       y1     |s|
|..............|.|..............|.|
~~~~~


To decompress a point take the value of y, and compute
(y^2-1)\*InvSqrt((y^2-1)\*(dy^2-1)), where InvSqrt(z)=1/sqrt(z).  This quantity
is one possible x value, its negation is the other. The top bit of the
compressed representation is 0 if the smaller one under the defined ordering
above is intended, and 1 otherwise. An algorithm for computing InvSqrt is
contained in the appendix.

This point compression format is from {{SchnorrQ}}, and the similar algorithm
there MAY be used instead to compute the x coordinates. Any method to decompress
points MAY be used provided it computes the correct answers. We call the
operation of compressing a point P into 32 bytes Compress(P), and decompression
Expand(S). Expand(Compress(P))=P for all P on the curve, and
Compress(Expand(S))=S if and only if S is a valid representation of a point.

Not all 32 byte strings represent valid points. Implementations MUST reject
invalid strings and check that decompression is successful. Strings are invalid
if they are not possible outputs of the compression operator.  In particular the
values of y0 and y1 MUST be less then p.

# Scalar multiplication

Below, we present two algorithms for scalar multiplication on the above curve:
one uses a simple fixed-window exponentiation without exploiting endomorphisms,
while the other uses endomorphisms to accelerate computation.  The execution of
operations in both algorithms has a regular pattern in order to enable
constant-time implementations and protect against timing and simple side
channel attacks.  Both algorithms use the same addition and doubling formulas.

First, we discuss explicit formulas and efficient projective coordinate
representations.

## Alternative Point Representations and Addition Laws

We use coordinates based on extended twisted Edwards coordinates introduced in
{{TwistedRevisited}}: the tuple (X, Y, Z, T) with Z nonzero and Z \* T = X \* Y
corresponds to a point (x, y) satisfying x = X/Z and y = Y/Z. The neutral point
in this representation is (0,1,1,0). The following slight variants are used in
the optimized scalar multiplication algorithm in order to save computations:
point representation R1 is given by (X,Y,Z,Ta,Tb), where T=Ta\*Tb;
representation R2 is (N, D, E, F) = (X+Y,Y-X,2Z,2dT); representation R3 is (N,
D, Z, T) = (X+Y,Y-X,Z,T); and representation R4 is (X,Y,Z). Similar "caching"
techniques were discussed in {{TwistedRevisited}} to accelerate repeated
additions of the same point. Converting between these representations is
straightforward.

A point doubling (DBL) takes an R4 point and produces an R1 point. For addition,
we first define an operation ADD_core that takes an R2 and R3 point and produces
an R1 point. This can be used to implement an operation ADD which takes an R1
and R2 point as inputs, converts the R1 point to R3, and then executes ADD_core.
Exposing these operations and the multiple representations helps save time by
avoiding redundant computations: the conversion of the first argument to ADD can
be done once if the argument will be used in multiple additions.

Below, we list the explicit formulas for the required point operations. These
formulas, which are adapted from {{Twisted}} and {{TwistedRevisited}}, are
complete: they have no exceptional cases, and therefore can be used in any
algorithm for computing scalar multiples without worrying about exceptional
procedure attacks {{Exceptional}}. Note that we do not explicitly note the point
format every time an addition or doubling is used, and assume that conversions
are done when required.

Doubling is computed as follows:

~~~~
DBL(X1, Y1, Z1):
  A = X1^2
  B = Y1^2
  C = 2*Z1^2
  D = A+B
  E = (X1+Y1)^2-D
  F = B-A
  G = C-F
  X3 = E*G
  Y3 = D*F
  Ta3 = E
  Tb3 = D
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

## Multiplication without endomorphisms

We begin by taking our input point P, and computing a table of points containing
T[0] = [1]P, T[1] = [3]P, ... , T[7] = [15]P as follows:

~~~~
Q = DBL(P)
Convert Q to R2 form
T[0] = P
Convert T[0] to R2 form
for i=1 to 7:
    T[i] = ADD_core(Q,T[i-1])
    Convert T[i] to R2 form
~~~~

Next, take m and reduce it modulo N.  Then, add N if necessary to ensure that m
is odd. At this point we recode m into a signed digit representation consisting
of 63 signed, odd digits d[i] in base 16. The following algorithm accomplishes
this task.

~~~~
for i=0 to 61:
    d[i] = (m mod 32) - 16
    m = (m - d[i]) / 16
d[62] = m
~~~~

Finally, the computation of the multiplication is as follows.

~~~~
Let ind = (abs(d[62])-1)/2
Let sign = sgn(d[62])
Q = sign*T[ind]
Convert Q into R4 form
for i from 61 to 0:
    Q = DBL(Q)
    Q = DBL(Q)
    Q = DBL(Q)
    Q = DBL(Q)
    ind = (abs(d[i])-1)/2
    sign = sgn(d[i])
    S = sign*T[ind]
    Q = ADD(Q, S)
return Q = (X/Z, Y/Z)
~~~~

As sign is either -1 or 1, the multiplication sign*T[ind] is simply a
conditional negation.  To negate a point (N, D, E, F) in R2 form one computes
(D, N, E, -F). The table lookups and conditional negations must be carefully
implemented as described in ``Security Considerations'' to avoid side-channel
attacks.  This algorithm MUST NOT be applied to points which are not N-torsion
points: it will produce the wrong answer.

## Multiplication with endomorphisms

This algorithm makes use of the identity [m]\*P = [a_1]\*P + [a_2]\*phi(P) +
[a_3]\*psi(P) + [a_4]\*psi(phi(P)), where a_1, a_2, a_3, and a_4 are 64-bit
scalars that depend on m (which, being less than N, has ~246 bits). The
multiexponentiation can be computed using a small table of 8 precomputed points
and 64 doublings and additions. This is considerably fewer operations than the
algorithm above, at the cost of a slightly more complicated implementation.

We describe each phase of the computation separately: the computation of the
endomorphisms, the scalar decomposition and recoding, the creation of the table
of precomputed points and, lastly, the computation of the final results. Each
section refers to constants listed in an appendix in order of appearance.

### Endomorphisms and Isogenies

The two endomorphisms phi and psi used to accelerate multiplication are computed
as phi(Q) = tau_dual(upsilon(tau(Q)) and psi(Q) = tau_dual(chi(tau(Q))).  Below,
we present procedures for tau, tau_dual, upsilon and chi, adapted from
{{FourQlib}}. Tau_dual produces an R1 point, while the other procedures produce
R4 points.

Nota Bene: tau produces points on a different curve, while upsilon and chi are
endomorphisms of that different curve. Tau and tau_dual are the isogenies
mentioned in the mathematical background above.  As a result the intermediate
results do not satisfy the equations of the curve E.  Implementers who wish to
check the correctness of these intermediate results are referred to {{Curve4Q}}.

~~~~
tau(X1, Y1, Z1):
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
  Ta2 = B-A
  D = 2*Z1^2-Ta2
  Tb2 = ctaudual1*X1*Y1
  X2 = Tb2*C
  Y2 = D*Ta2
  Z2 = D*C
return(X2, Y2, Z2, Ta2, Tb2)
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

### Table Precomputation

This stage consists in computing a table of 8 points in representation R2 as
shown below.  Computations Q = psi(P), R = phi(P) and S = psi(phi(P)) are
carried out using formulas from the previous section.

~~~~~
Q is phi(P)
R is psi(P)
S is psi(Q)
Convert Q to R3
Convert R to R3
Convert S to R3
T[0] is P in R2
T[1] is ADD_core(Q, T[0])  (P+Q)
Convert T[1] to R2
T[2] is ADD_Core(R, T[0])  (P+R)
Convert T[2] to R2
T[3] is ADD_Core(R, T[1])  (P+Q+R)
Convert T[3] to R2
T[4] is ADD_Core(S, T[0])  (P+S)
Convert T[4] to R2
T[5] is ADD_Core(S, T[1])  (P+Q+S)
Convert T[5] to R2
T[6] is ADD_Core(S, T[2])  (P+R+S)
Convert T[6] to R2
T[7] is ADD_Core(S, T[3])  (P+Q+R+S)
Convert T[7] to R2
~~~~~

### Scalar Decomposition and Recoding

This stage has two parts. The first one consists in decomposing the scalar into
four 64-bit integers, and the second one consists in recoding these integers
into a form that can be used to efficiently and securely compute the scalar
multiplication.

The decomposition step uses four fixed vectors called b1, b2, b3, b4, with 64
bit entries each.  In addition, we have l1, l2, l3, l4, which are 192-bit
integers used to implement rounding.  All these values are listed in the
appendix.

First, compute c = 2\*b1 - b2 + 5\*b3 + 2\*b4 and c' = 2\*b1 - b2 + 5\*b3 + b4.
Note that c and c' can be precalculated.  Next, compute ti = floor(li\*m/2^256)
for i between 1 and 4, and then compute a = (a1, a2, a3, a4) = (m,0,0,0) -
t1\*b1 - t2\*b2 - t3\*b3 - t4\*b4.  Precisely one of a+c and a+c' has an odd
first coordinate: this is the vector v that is fed into the scalar recoding
step. Note that the entries of this vector are 64 bits, so intermediate values
in the calculation above can be truncated to this width.

The recoding step takes the vector v=(v1, v2, v3, v4) from the previous step and
outputs two arrays m[0]..m[64] and d[0]..d[64]. Each entry of d is between 0 and
7, and each entry in m is -1 or 0. The recoding algorithm is detailed below.
bit(x, n) denotes the nth bit of x, counting from least significant to most,
starting with 1.

~~~~~
m[64]=-1
for i=0 to 63 do:
   d[i] = 0
   m[i] = -bit(v1, i+1)
   for j = 2 to 4 do:
      d[i] = d[i]+bit(vj, 0)*2^(j-2)
      c = (bit(v1, i+1) or bit(vj,0)) xor bit(v1, i+1)
      vj = vj/2+c
d[64] = v2+2*v3+4*v4
~~~~~

### Final Computation

We now describe the last step in the endomorphism based algorithm for computing
scalar multiplication. On inputs m and P, the algorithm first computes the
precomputed table T with 8 points (see {{table-precomputation}}) and then carries
out the scalar decomposition and scalar recoding to produce the two arrays
m[0]..m[64] and d[0]..d[64].  See {{scalar-decomposition-and-recoding}}.

Define s[i] to be 1 if m[i] is -1 and -1 if m[i] is 0. Then the multiplication
is completed by the following pseudocode:

~~~~~
Q = s[64]*T[d[64]]
Convert Q to R4
for i=63 to 0 do:
    Q = DBL(Q)
    Q = ADD(Q, s[i]*T[di])
return Q = (X/Z, Y/Z)
~~~~~

Multiplication by s[i] is simply a conditional negation. To negate an R2 point
(N, D, E, F) one computes (D, N, E ,-F). It is important to do this (as well as
the table lookup) in constant time, i.e., the execution of branches and memory
accesses MUST NOT depend on secret values (see ``Security Considerations'' for
more details).

The optimized multiplication algorithm above only works properly for N-torsion
points. Implementations MUST NOT use this algorithm on anything that is not
known to be an N-torsion point. Otherwise, it will produce the wrong answer,
with extremely negative consequences for security.

# Use of the scalar multiplication primitive for Diffie-Hellman Key Agreement

The above scalar multiplication algorithms can be used to implement
Diffie-Hellman with cofactor.

~~~~
DH(m, P):
      Ensure P on curve and if not return FAILURE
      Q = [392]*P
      Compute [m]*Q
Return [m]*Q in affine coordinates
~~~~

The role of the separate multiplication by 392 is to ensure that Q is an
N-torsion point so that the scalar multiplication algorithms above may be used
safely to produce correct results.  In other words, as the cofactor is greater
than one, Diffie-Hellman computations using Curve4Q MUST always use cofactor
clearing (as defined above). Note that the multiplication by the cofactor 392
can be safely computed with nine doublings and two additions via the double and
add method.

[[ RLB: Why not just transform the coefficient instead? ]]

Two users, Alice and Bob, can carry out the following steps to derive a shared
key: each picks a random string of 32 bytes, mA and mB, respectively. Alice
computes the public key A = Compress(DH(mA, G)), and Bob computes the public key
B = Compress(DH(mB, G)). They exchange A and B, and then Alice computes KAB =
DH(mA, Expand(B)) while Bob computes KBA = DH(mB, Expand(A)), which produces the
shared point K = KAB = KBA. The y coordinate of K, represented as a 32 byte
string as detailed in {{representation-of-curve-points}} is the shared
secret.

If the received strings are not valid points, the DH function has failed to
compute an answer. Implementations SHOULD return a random 32 byte string as well
as return an error, to prevent bugs when applications ignore return codes. They
MUST signal an error when decompression fails.

Implementations MAY use any method to carry out these calculations, provided
that it agrees with the above function on all inputs and failure cases, and does
not leak information about secret keys. For example, refer to the constant-time
fixed-base scalar multiplication algorithm implemented in {{FourQlib}} to
accelerate the computation of DH(m, G).

# IANA Considerations

[RFC Editor: please remove this section prior to publication]
This document has no IANA actions.

# Security Considerations

The best known algorithms for the computation of discrete logarithms on Curve4Q
are parallel versions of the Pollard rho algorithm in [Distinguished]. On
Curve4Q these attacks take on the order of 2^123 group operations to compute a
single discrete logarithm. The additional endomorphisms have large order, and so
cannot be used to accelerate generic attacks. Quadratic fields are not affected
by any of the index calculus attacks used over larger extension fields.

Implementations MUST check that input points properly decompress to points on
the curve. Removing such checks may result in extremely effective attacks. The
curve is not twist-secure: implementations using single coordinate ladders MUST
validate points before operating on them. In the case of protocols that require
contributory behavior, when the identity is the output of the DH primitive it
MUST be rejected and failure signaled to higher levels. Notoriously {{?RFC5246}}
without {{?RFC7627}} is such a protocol.

Implementations MUST ensure that execution of branches and memory addresses
accessed do not depend on secret data.  The time variability introduced by
secret-dependent operations have been exploited in the past via timing and cache
attacks to break implementations.  Side-channel analysis is a constantly moving
field, and implementers must be extremely careful to ensure that operations do
not leak any secret information. Using ephemeral private scalars for each
operation (ideally, limiting the use of each private scalar to one single
operation) can reduce the impact of side-channel attacks.  However, this might
not be possible for many applications of Diffie-Hellman key agreement.

In the future quantum computers may render the discrete logarithm problem easy
on all abelian groups through Shor's algorithm. Data intended to remain
confidential for significantly extended periods of time SHOULD NOT be protected
with any primitive based on the hardness of factoring or the discrete log
problem (elliptic curve or finite field).

--- back

# Constants

ctau1 = 221360928884514619410\*i + 33754435779700894835198039471158097091

ctaudual1 = 170141183460469231510326374831369486353\*i+ 99231301967130569661901792840482943028

cphi0 = 49615650983565284830950896420241471514\*i + 110680464442257309687

cphi1 = 131306912742858181648727312260439119609\*i + 92233720368547758087

cphi2 = 160666015865631300014011952927357137809\*i + 276701161105643274261

cphi3 = 107027644557995218531204623577807990436\*i + 36893488147419103235

cphi4 = 24279268184862963117522688682631129173\*i + 55340232221128654851

cphi5 = 92472642025247131565767320804994133491\*i + 184467440737095516175

cphi6 = 14804100590025031399847337894104161255\*i + 332041393326771929112

cphi7 = 76283848507754718862858058709728786458\*i + 442721857769029238819

cphi8 = 41635071732389019719735756359456329456\*i + 3135946492530623774960

cphi9 = 21045324596686230484035983431638590725\*i + 39844967199212631493615

cpsi1 = 4095177184363520459066\*i + 57123674603396429897431647433607300847

cpsi2 = 44824135016688633386011024159913800562\*i + 4205857648805777768771

cpsi3 = 101947809620085063283442671593521101409\*i + 110680464442257309705

cpsi4 = 68193373840384168448244632122363004318\*i + 170141183460469231621006839273626796022

d = 0x5e472f846657e0fcb3821488f1fc0c8d\*i + 0x00000000000000e40000000000000142

l1 = 50127518246259276682880317011538934615153226543083896339791

l2 = 22358026531042503310338016640572204942053343837521088510715

l3 = 5105580562119000402467322500999592531749084507000101675068

l4 = 19494034873545274265741574254707851381713530791194721254848

b1 = [650487742939046294, -1397215820276968864, 523086274270593807,
   -598824378691085905]

b2 = [2110318963211420372, -1, 1, 2727991412926801872]

b3 = [1705647224544756482, 199320682881407569,
   -3336360048424633503, 765171327772315031]

b4 = [1400113754146392127, 3540637644719456050,
-471270406870313397, -1789345740969872106]

Gx = 0x1E1F553F2878AA9C96869FB360AC77F6\*i + 0x1A3472237C2FB305286592AD7B3833AA

Gy = 0x6E1C4AF8630E024249A7C344844C8B5C\*i + 0x0E3FEE9BA120785AB924A2462BCBB287

# Inversion and Square roots

Inversion of nonzero elements of GF(p) can be computed in constant-time using
one exponentiation via Fermat's Little Theorem: 1/a = a^(p - 2) = a^(2^127 - 3).

The following algorithm for computing inverse square roots in GF(p^2) is an
adaptation of Algorithm 8 from {{SQRT}}. Note that (p-3)/4 is 2^125-1, and there
is a very short addition chain to compute this value.

~~~~~
   InvSqrt(a+b*i):
        if b = 0:
            t = a^((p-3)/4)
            if a*t^2=1:
               return t+0*i
            else:
                return 0 + t*i
        else:
            n = a^2+b^2
            s = (n)^((p-3)/4)
            if n*s^2 is not 1:
                return FAILURE
            c = n*s
            delta = (a+c)/2
            g = (delta)^((p-3)/4)
            if delta*g^2 = 1:
                x0 = delta*g
                x1 = (b/2)*delta
             else:
                 x1 = delta*g
                 x0 = (b/2)*delta
             return (x0+i*x1)*s
~~~~~
