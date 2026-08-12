"""Faithful Python port of the Numerical Recipes eigensolver used by ``qcdf.f``.

Why this exists
---------------
``qcdf.f`` adds the free part of the Hamiltonian to the *diagonal only* of the
orthonormal-basis matrix::

    HNU(I1,I1) = (1 - lambda^2) * rmq^2 * HNU0(I1) + HNU(I1,I1)

That step is **basis dependent**.  ``Z^T H0 Z`` is not diagonal (its largest
off-diagonal element is ~0.8 at N=3, B=1, 2K=21), so discarding the off-diagonal
part makes the resulting spectrum depend on which eigenvectors the norm-matrix
diagonalization happened to return.  The norm matrix is highly degenerate -- 26
distinct eigenvalues among 189 states -- so that choice is wide open, and LAPACK
and Numerical Recipes make different ones.  Feeding the two different ``Z``
matrices through the same code changes M^2 by up to 5.8.

Consequences, measured at N=3, NF=1, B=1, 2K=21, lambda=0.3325:

======================================  ==========================
assembly                                spectrum
======================================  ==========================
full ``Z^T H0 Z``                       basis independent to 6e-13
diagonal of ``Z^T H0 Z`` (what qcdf.f   basis dependent, differs by
does)                                   up to 5.8 between Z choices
======================================  ==========================

The physically correct ground state is 10.390824; ``qcdf.f`` reports 10.390380,
a 4.3e-5 relative shift.  That is far below the precision the paper quotes
(Table I gives M_bar/g = 10.71(2), i.e. 0.2%), so the published physics stands.
But reproducing the historical numbers *exactly* requires reproducing the
historical eigenvectors, which is what this module provides.

``dlcq.read_python`` exposes both: ``assembly="fortran"`` for bit-level
regression against the 1990 output, ``assembly="exact"`` for new work.

The three routines mirror ``TRR8`` (Numerical Recipes ``tred2``), ``TQR8``
(``tqli``) and ``ESRTR8`` (an ascending selection sort) line for line, including
the ``ABS(E(M))+DD .EQ. DD`` convergence test whose exact floating-point
behaviour determines the result.
"""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
    _HAVE_NUMBA = True
except ImportError:                                    # pragma: no cover
    _HAVE_NUMBA = False

    def njit(*args, **kwargs):
        def wrap(fn):
            return fn
        return wrap(args[0]) if args and callable(args[0]) else wrap


__all__ = ["nr_eigh", "tred2", "tqli", "eigsrt"]


@njit(cache=True)
def _sign(a, b):
    """Fortran ``SIGN(a, b)``: magnitude of ``a``, sign of ``b`` (``b == 0`` -> +)."""
    return abs(a) if b >= 0.0 else -abs(a)


@njit(cache=True)
def tred2(a, n):
    """Householder reduction to tridiagonal form (Numerical Recipes ``tred2``).

    ``a`` is a 1-based ``(n+1, n+1)`` working array, modified in place to hold
    the orthogonal transformation.  Returns ``(d, e)``, also 1-based.
    """
    d = np.zeros(n + 1)
    e = np.zeros(n + 1)

    if n > 1:
        for i in range(n, 1, -1):
            l = i - 1
            h = 0.0
            scale = 0.0
            if l > 1:
                for k in range(1, l + 1):
                    scale += abs(a[i, k])
                if scale == 0.0:
                    e[i] = a[i, l]
                else:
                    for k in range(1, l + 1):
                        a[i, k] /= scale
                        h += a[i, k] * a[i, k]
                    f = a[i, l]
                    g = -_sign(np.sqrt(h), f)
                    e[i] = scale * g
                    h -= f * g
                    a[i, l] = f - g
                    f = 0.0
                    for j in range(1, l + 1):
                        a[j, i] = a[i, j] / h
                        g = 0.0
                        for k in range(1, j + 1):
                            g += a[j, k] * a[i, k]
                        if l > j:
                            for k in range(j + 1, l + 1):
                                g += a[k, j] * a[i, k]
                        e[j] = g / h
                        f += e[j] * a[i, j]
                    hh = f / (h + h)
                    for j in range(1, l + 1):
                        f = a[i, j]
                        g = e[j] - hh * f
                        e[j] = g
                        for k in range(1, j + 1):
                            a[j, k] -= f * e[k] + g * a[i, k]
            else:
                e[i] = a[i, l]
            d[i] = h

    d[1] = 0.0
    e[1] = 0.0
    for i in range(1, n + 1):
        l = i - 1
        if d[i] != 0.0:
            for j in range(1, l + 1):
                g = 0.0
                for k in range(1, l + 1):
                    g += a[i, k] * a[k, j]
                for k in range(1, l + 1):
                    a[k, j] -= g * a[k, i]
        d[i] = a[i, i]
        a[i, i] = 1.0
        if l >= 1:
            for j in range(1, l + 1):
                a[i, j] = 0.0
                a[j, i] = 0.0
    return d, e


@njit(cache=True)
def tqli(d, e, n, z):
    """QL with implicit shifts on a tridiagonal matrix (Numerical Recipes ``tqli``).

    ``d``, ``e`` are 1-based; ``z`` is the 1-based transformation from
    :func:`tred2`, updated in place to hold the eigenvectors as columns.
    """
    if n > 1:
        for i in range(2, n + 1):
            e[i - 1] = e[i]
        e[n] = 0.0

        for l in range(1, n + 1):
            iter_count = 0
            while True:
                m = n
                for mm in range(l, n):
                    dd = abs(d[mm]) + abs(d[mm + 1])
                    # The exact NR convergence test; do not "simplify" it.
                    if abs(e[mm]) + dd == dd:
                        m = mm
                        break
                if m == l:
                    break
                if iter_count == 1000:
                    raise ValueError("tqli: too many iterations")
                iter_count += 1
                g = (d[l + 1] - d[l]) / (2.0 * e[l])
                r = np.sqrt(g * g + 1.0)
                g = d[m] - d[l] + e[l] / (g + _sign(r, g))
                s = 1.0
                c = 1.0
                p = 0.0
                for i in range(m - 1, l - 1, -1):
                    f = s * e[i]
                    b = c * e[i]
                    if abs(f) >= abs(g):
                        c = g / f
                        r = np.sqrt(c * c + 1.0)
                        e[i + 1] = f * r
                        s = 1.0 / r
                        c *= s
                    else:
                        s = f / g
                        r = np.sqrt(s * s + 1.0)
                        e[i + 1] = g * r
                        c = 1.0 / r
                        s *= c
                    g = d[i + 1] - p
                    r = (d[i] - g) * s + 2.0 * c * b
                    p = s * r
                    d[i + 1] = g + p
                    g = c * r - b
                    for k in range(1, n + 1):
                        f = z[k, i + 1]
                        z[k, i + 1] = s * z[k, i] + c * f
                        z[k, i] = c * z[k, i] - s * f
                d[l] -= p
                e[l] = g
                e[m] = 0.0
    return d, z


@njit(cache=True)
def eigsrt(d, v, n):
    """Ascending selection sort of eigenvalues with matching columns (``ESRTR8``)."""
    for i in range(1, n):
        k = i
        p = d[i]
        for j in range(i + 1, n + 1):
            if d[j] < p:
                k = j
                p = d[j]
        if k != i:
            d[k] = d[i]
            d[i] = p
            for j in range(1, n + 1):
                tmp = v[j, i]
                v[j, i] = v[j, k]
                v[j, k] = tmp
    return d, v


def nr_eigh(matrix):
    """Symmetric eigendecomposition exactly as ``qcdf.f``'s ``DIAG`` performs it.

    Drop-in for ``scipy.linalg.eigh`` -- returns ``(w, z)`` with ascending ``w``
    and eigenvectors in the columns of ``z`` -- but reproduces the Numerical
    Recipes result, including its particular choice of basis inside degenerate
    eigenspaces.  That choice is what makes ``qcdf.f``'s diagonal-only free-part
    assembly reproducible.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    n = matrix.shape[0]
    if n == 0:
        return np.zeros(0), np.zeros((0, 0))

    a = np.zeros((n + 1, n + 1), dtype=np.float64)
    a[1:, 1:] = matrix

    d, e = tred2(a, n)
    d, a = tqli(d, e, n, a)
    d, a = eigsrt(d, a, n)

    return d[1:].copy(), a[1:, 1:].copy()
