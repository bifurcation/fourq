#!/usr/bin/python

from random import getrandbits
from time import time

# GF(p^2) arithmetic

p1271 = (1 << 127) - 1

def fp1271add(a, b):
    return (a + b) % p1271

def fp1271sub(a, b):
    return (a - b) % p1271

def fp1271mul(a, b):
    return (a * b) % p1271

def fp2_1271mul(a, b):
    a0b0 = fp1271mul(a[0], b[0])
    a1b0 = fp1271mul(a[1], b[0])
    a0b1 = fp1271mul(a[0], b[1])
    a1b1 = fp1271mul(a[1], b[1])
    return (fp1271sub(a0b0, a1b1), fp1271add(a0b1, a1b0))


# GF(p)

p25519 = (1 << 255) - 19

def fp25519mul(x, y):
    return (x * y) % p25519


# Speed test

TEST_SAMPLES = 1000
print "Samples: {}". format(TEST_SAMPLES)

corpus = [getrandbits(256) for i in range(TEST_SAMPLES)]

tic = time()
corpus1271 = [(x & p1271, (x >> 128) % p1271) for x in corpus]
out1271 = range(TEST_SAMPLES)
for i in range(0, TEST_SAMPLES-1):
    out1271[i] = fp2_1271mul(corpus1271[i], corpus1271[i+1])
toc = time()
print "GF(p_1271^2) = {:3.2f}ms".format((toc - tic)*1000)

tic = time()
corpus25519 = [x % p25519 for x in corpus]
out25519 = range(TEST_SAMPLES)
for i in range(0, TEST_SAMPLES-1):
    out25519[i] = fp25519mul(corpus25519[i], corpus25519[i+1])
toc = time()
print "GF(p_25519) = {:3.2f}ms".format((toc - tic)*1000)
