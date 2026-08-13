# Multiple flavours (NF > 1)

Flavour was added to `qcdf.f` on 5/8/90 — the second-to-last dated change in its
header — and **no figure in the paper uses it**. The path therefore shipped with
no coverage, and this document is what checking it turned up.

It matters because it is the obvious direction for new work: unequal quark masses
give the analogue of a kaon, and of a strange baryon, neither of which the paper
computes.

## Driving it

Both solvers take one mass per flavour and one flavour quantum number per
flavour. `qcdf.f` fixes `rmq(1) = 1` and reads `rmq(2..NF)` as **ratios** to it;
`iflv(f)` is the **net quark number of flavour f** (quarks minus antiquarks), so
a meson is flavour-neutral overall and a baryon's numbers must sum to `N·B`.

```bash
# flavour-neutral meson, second flavour twice as heavy
bash fortran/run_case.sh runs/nf2 3 2 0 0.3325 -1.0 0 12   2.0   0 0
```

```python
from dlcq.read_python import run_python
r = run_python(N=3, NF=2, B=0, K_code=12, rlamb=0.3325,
               rmq=[1.0, 2.0], iflv=[0, 0])       # neutral
r = run_python(N=3, NF=2, B=1, K_code=13, rlamb=0.3325,
               rmq=[1.0, 2.0], iflv=[2, 1])       # "uus" baryon
```

`FortranProvider.get` takes the same `rmq`/`iflv`. `NF ≤ 3`: `qcdf.f` dimensions
`rmq(03)` and `iflv(03)`.

Note B = 1, N = 3 needs an **odd** `K_code` — three quarks of odd momentum — and
asking for an even one returns zero states rather than an error.

## Defect: `qcdf.f` cannot reach a flavour non-singlet

Line 253, immediately before `CALL QCDSTA` and carrying the author's own
`??????` marker:

```fortran
      iflv(1) = N*B
      CALL QCDSTA
```

is **unconditional**. Whatever was entered at the `input flavor qm number for
flavor 1` prompt is discarded. Only the channel with `iflv(1) = N·B` is
reachable; any other asks for a total quark number the code cannot satisfy and
generates **zero states rather than an error**.

So the Fortran cannot compute a flavour non-singlet — no kaon analogue, and no
flavoured baryon except the one whose first flavour carries all the baryon
number. Asking for `iflv = [1, −1]` at 2K = 12 gives 0 states where the Python
solver gives 316.

`dlcq/read_python.py` honours the value it is given, and only defaults
`iflv[0] = N·B` when the caller supplies nothing. `fortran/qcdf.f` is left
unpatched, as everywhere else in this repository.

## The colour-array overflow arrives much sooner

`qcdf.f` needs `2L + 4` colour-index slots for a matrix element between
`L`-parton states and dimensions that array with 25, so runs reaching 11 partons
are silently corrupted (`fortran-color-overflow.md`). Pauli exclusion allows N
quarks per *(momentum, flavour)*, so a second flavour packs more partons into the
same momentum — and the overflow arrives at half the K:

| NF | first overflowing 2K | max partons there |
|---|---|---|
| 1 | 24 | 12 |
| **2** | **12** | **12** |

Anything past 2K = 10 at NF = 2 therefore needs either the widened build
(`fortran/widen.py`, which raises the caps to 54) or the Python solver, which is
unaffected. This is a much tighter constraint than at NF = 1 and is easy to miss:
the symptom is non-positive M² at the *bottom* of the spectrum, exactly where a
ground-state mass is read.

## Validation

There is nothing published to compare against, so the checks are identities,
inequalities, and cross-solver agreement. All are in `tests/test_flavour.py`.

**Exact identities.** Both sum rules hold to 10⁻¹² across seven configurations —
NF = 1 and 2, mesons and baryons, degenerate and unequal masses, neutral and
non-singlet. The momentum rule is the sharper of the two: every basis state
carries total momentum K regardless of how its partons are flavoured, so it holds
state by state and directly tests that the enlarged basis is built and weighted
correctly.

**Flavour symmetry.** At `m₁ = m₂` the neutral channel (`iflv = [0,0]`, 443
states) and the non-singlet channel (`[1,−1]`, 316 states) must be degenerate.
They are, to 10⁻⁹ — built from different numbers of basis states, so this tests
the flavour bookkeeping rather than restating it.

**Variational inequalities.** NF = 2 contains every NF = 1 configuration, so its
ground state must sit at or below the NF = 1 one, even with a heavier second
quark. It does. And doubling the mass of a quark the channel actually contains
must raise the state — checked for the non-singlet meson and the *uus* baryon.

**Both solvers.** Where `qcdf.f` can run — the neutral channel — the two agree
at 2K = 12 on state generation (510 → 443, identically) and on the ground state
to **3 × 10⁻⁷**, the same standard as the NF = 1 path and far inside the ~10⁻⁴
floor of `basis-dependence.md`.

## What becomes computable

| state | `iflv` | 2K | M² |
|---|---|---|---|
| neutral meson, degenerate | [0, 0] | 12 | 4.468521 |
| neutral meson, m₂ = 2m₁ | [0, 0] | 12 | 4.472366 |
| non-singlet meson, degenerate | [1, −1] | 12 | 4.468521 |
| **non-singlet meson, m₂ = 2m₁** | [1, −1] | 12 | **9.272800** |
| **baryon *uud*** | [2, 1] | 13 | **10.166076** |
| **baryon *uus*, m₂ = 2m₁** | [2, 1] | 13 | **16.821694** |

The bolded rows are physics the paper does not contain. They are small-K numbers
quoted to show the machinery runs, not results — a real calculation would sweep K
and extrapolate as Table I does, which `_extrapolated_mass` already supports.
