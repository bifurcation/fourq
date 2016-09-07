#!/usr/bin/env python

from random import getrandbits
from time import time

from fields import GFp2, GFp25519, p1271, p25519
import curve4q
import curve25519

# Adjust these if you want more/fewer samples
FIELD_TEST_LOOPS = 1000
DH_TEST_LOOPS = 100

def compare_fields():
    base_corpus = [getrandbits(256) for i in range(FIELD_TEST_LOOPS)]
    corpus1271 = [(x % p1271, (x >> 128) % p1271) for x in base_corpus]
    corpus25519 = [x % p25519 for x in base_corpus]

    def speedtest(f):
        # Test in GFp2
        tic = time()
        for i in range(0, FIELD_TEST_LOOPS-1):
            f(GFp2, i, corpus1271)
        toc = time()
        t1271 = 1000 * (toc - tic)

        # Test in GFp25519
        tic = time()
        for i in range(0, FIELD_TEST_LOOPS-1):
            f(GFp25519, i, corpus25519)
        toc = time()
        t25519 = 1000 * (toc - tic)

        return (t1271, t25519)

    tests = {
        "add": lambda field, i, corpus: field.add(corpus[i], corpus[i+1]),
        "mul": lambda field, i, corpus: field.mul(corpus[i], corpus[i+1]),
        "sqr": lambda field, i, corpus: field.sqr(corpus[i]),
        "inv": lambda field, i, corpus: field.inv(corpus[i])
    }

    print "===== Time for {} field operations =====".format(FIELD_TEST_LOOPS)
    print
    print "{:5s}   {:>8s}   {:>8s}".format("Op", "GFp2", "GFp25519")
    for name in tests:
        (t1271, t25519) = speedtest(tests[name])
        print "{:5s} {:8.2f}ms {:8.2f}ms".format(name, t1271, t25519)
    print

def compare_ops():
    opcounts = {}

    m = getrandbits(256)
    G = curve4q.AffineToR1(curve4q.Gx, curve4q.Gy)
    k = '77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a'.decode('hex')
    u = '0900000000000000000000000000000000000000000000000000000000000000'.decode('hex')

    G392 = curve4q.MUL_endo(392, G)
    T_windowed = curve4q.table_windowed(G)
    T_endo = curve4q.table_endo(G)
    T392_windowed = curve4q.table_windowed(G392)
    T392_endo = curve4q.table_endo(G392)

    # R1toR2
    GFp2.ctr_reset()
    Q = curve4q.R1toR2(G)
    opcounts["R1toR2"] = GFp2.ctr()

    # R1toR3
    GFp2.ctr_reset()
    Q = curve4q.R1toR3(G)
    opcounts["R1toR3"] = GFp2.ctr()

    # R2toR4
    G2 = curve4q.R1toR2(G)
    GFp2.ctr_reset()
    Q = curve4q.R2toR4(G2)
    opcounts["R2toR4"] = GFp2.ctr()

    # ADD_core
    P = curve4q.R1toR3(G)
    Q = curve4q.R1toR2(G)
    GFp2.ctr_reset()
    R = curve4q.ADD_core(P, Q)
    opcounts["ADD_core"] = GFp2.ctr()

    # ADD
    P = G
    Q = curve4q.R1toR2(G)
    GFp2.ctr_reset()
    R = curve4q.ADD(P, Q)
    opcounts["ADD"] = GFp2.ctr()

    # DBL
    GFp2.ctr_reset()
    Q = curve4q.DBL(G)
    opcounts["DBL"] = GFp2.ctr()

    # MUL_windowed
    GFp2.ctr_reset()
    Q = curve4q.MUL_windowed(m, G)
    opcounts["MUL_windowed"] = GFp2.ctr()

    # MUL_windowed_fixed
    GFp2.ctr_reset()
    Q = curve4q.MUL_windowed(m, G, table=T_windowed)
    opcounts["MUL_windowed_fixed"] = GFp2.ctr()

    # Phi
    GFp2.ctr_reset()
    phiP = curve4q.phi(G)
    opcounts["phi"] = GFp2.ctr()

    # Psi
    GFp2.ctr_reset()
    psiP = curve4q.psi(G)
    opcounts["psi"] = GFp2.ctr()

    # MUL_endo
    GFp2.ctr_reset()
    mP = curve4q.MUL_endo(m, G)
    opcounts["MUL_endo"] = GFp2.ctr()

    # MUL_endo_fixed
    GFp2.ctr_reset()
    Q = curve4q.MUL_endo(m, G, table=T_endo)
    opcounts["MUL_endo_fixed"] = GFp2.ctr()

    # DH_windowed
    GFp2.ctr_reset()
    mP = curve4q.DH_windowed(m, G[:2])
    opcounts["DH_windowed"] = GFp2.ctr()

    # DH_windowed_fixed
    GFp2.ctr_reset()
    mP = curve4q.DH_windowed(m, G[:2], table=T392_windowed)
    opcounts["DH_windowed_fixed"] = GFp2.ctr()

    # DH_endo
    GFp2.ctr_reset()
    mP = curve4q.DH_endo(m, G[:2])
    opcounts["DH_endo"] = GFp2.ctr()

    # DH_endo_fixed
    GFp2.ctr_reset()
    mP = curve4q.DH_endo(m, G[:2], table=T392_endo)
    opcounts["DH_endo_fixed"] = GFp2.ctr()

    # x25519
    GFp25519.ctr_reset()
    ku = curve25519.x25519(k, u)
    opcounts["x25519"] = GFp25519.ctr()

    rows = ["R1toR2", "R1toR3", "R2toR4", "ADD_core", "ADD", "DBL",
            "phi", "psi", "psiphi",
            "MUL_windowed", "MUL_windowed_fixed", "MUL_endo", "MUL_endo_fixed",
            "DH_windowed", "DH_windowed_fixed", "DH_endo", "DH_endo_fixed",
            "x25519"]

    print "===== Field operation count ====="
    print
    print "{:18s} {:>7s} {:>7s} {:>7s} {:>7s}".format("", "M", "S", "A", "I")
    for name in rows:
        if name not in opcounts:
            continue
        opctr = opcounts[name]
        print "{:20s} {:7.1f} {:7.1f} {:7.1f} {:7.1f}".format(name, opctr[2], opctr[1], opctr[0], opctr[3])
    print

def compare_time():
    G = (curve4q.Gx, curve4q.Gy)
    u = '09'.ljust(64, '0').decode('hex')

    G392 = curve4q.MUL_endo(392, curve4q.AffineToR1(curve4q.Gx, curve4q.Gy))
    T392_windowed = curve4q.table_windowed(G392)
    T392_endo = curve4q.table_endo(G392)

    coeff4q = [getrandbits(256) for i in range(DH_TEST_LOOPS)]
    coeff25519 = ["{:x}".format(m).ljust(64, '0').decode('hex') for m in coeff4q]

    time4q_windowed = 0
    time4q_windowed_fixed = 0
    time4q_endo = 0
    time4q_endo_fixed = 0
    time25519 = 0
    for i in range(len(coeff4q)):
        tic = time()
        mP = curve4q.DH_windowed(coeff4q[i], G)
        toc = time()
        time4q_windowed += toc - tic

        tic = time()
        mP = curve4q.DH_windowed(coeff4q[i], G, table=T392_windowed)
        toc = time()
        time4q_windowed_fixed += toc - tic

        tic = time()
        mP = curve4q.DH_endo(coeff4q[i], G)
        toc = time()
        time4q_endo += toc - tic

        tic = time()
        mP = curve4q.DH_endo(coeff4q[i], G, table=T392_endo)
        toc = time()
        time4q_endo_fixed += toc - tic

        tic = time()
        ku = curve25519.x25519(coeff25519[i], u)
        toc = time()
        time25519 += toc - tic

    print "===== Time for {} field operations =====".format(DH_TEST_LOOPS)
    print
    print "{:<28s} {:>7.2f}ms".format("Curve4Q (windowed)", 1000 * time4q_windowed)
    print "{:<28s} {:>7.2f}ms".format("Curve4Q (win fixed base)", 1000 * time4q_windowed_fixed)
    print "{:<28s} {:>7.2f}ms".format("Curve4Q (endomorphisms)", 1000 * time4q_endo)
    print "{:<28s} {:>7.2f}ms".format("Curve4Q (endo fixed base)", 1000 * time4q_endo_fixed)
    print "{:<28s} {:>7.2f}ms".format("Curve25519", 1000 * time25519)

if __name__ == "__main__":
    compare_fields()
    compare_ops()
    compare_time()
