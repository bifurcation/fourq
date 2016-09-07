curve4q.py
==========

This directory contains an implementation of Curve4Q in python, following the
instructions in the specification as directly as possible.  You should be able
to map directly from the specification to the functions in `curve4q.py`.

For purposes of comparison, there is also an implementation of Curve25519.
Again, this implementation is derived as directly as possible from the
specification, namely [RFC 7748](https://tools.ietf.org/html/rfc7748).

Both implementations come with tests to verify their correctness.  (For Curve4Q,
this also helps verify the correctnes of the spec.)  The tests for Curve25519
come directly from the RFC.  The test for Curve4Q are taken from the Microsoft
[FourQlib](http://research.microsoft.com/en-us/projects/fourqlib/)
implementation.

You can get an idea of how the curves compare by running the included comparison
script:

```
> python curve4q.py     # Correctness tests
> python curve25519.py  # Correctness tests
> python compare.py     # Comparison tests
```

For the most part, these comparisons are pretty rough.  Using python makes the
code easier to read, but it means that we can't really optimize the field
arithmetic.  So the overall curve performance numbers could vary depending on
how the underlying field operations are implemented on a given platform.

However, the operation count numbers only really depend on the specification, so
they should be pretty reliable.

