# Independent periodic-image fixture certificates

`cases.json` is the fixed numerical input, in Angstrom. Cell rows are the
vectors a,b,c (columns of H), not inverse rows. Every fixture declares cube=16;
the oracle always enumerates every ordered pair and all x,y,z in [-16,16].
It filters disabled lattice axes and zero-self/self policy before the direct
scalar displacement and strict `(dx*dx+dy*dy)+dz*dz < R*R` predicate. It does
not import production geometry, use minimum images, or derive its cube from
production bounds. Sorting and comparing full (i,j,Sx,Sy,Sz) tuples is exact.

The 30 explicit snapshots include all eight masks of multi4 and mace3.
`self_disabled` and `self_short` cover both policy/radius variants. Boundary
positions use the immediate representable f64 neighbors of 2; they are never
sent to ASE. Empty input still goes through production validation in the
separate graph input tests.

## Completeness independent of the production search

For the triangular cells a=(A,0,0), b=(B,C,0), c=(D,E,F), the direct component
equations give successively

```
|Sz| <= (Q+max|delta_z|)/|F|
|Sy| <= (Q+max|delta_y|+|E|*bound_z)/|C|
|Sx| <= (Q+max|delta_x|+|B|*bound_y+|D|*bound_z)/|A|
Q = R + 1e-8
```

The +1e-8 is search-only: it never changes physical membership. The numeric
bounds below are rounded outward (bounds below an integer admit its floor).
They also exclude any missed image within the official 1e-8 cutoff margin.
These fixture constants do not involve an inverse or production certification.

```
unreduced: |Sz| < (4.8+1e-8)/40 < .120000001 => Sz=0;
|Sy| < (4.8+1e-8+.1)/10 < .490000002 => Sy=0;
|Sx| < (4.8+1e-8+19.7)/10 < 2.450000002 => |Sx|<=2.
oblique: Sz=0;
|Sy| < (4.8+1e-8+3)/10 < .780000002 => Sy=0;
|Sx| < (4.8+1e-8+6.37)/10 < 1.117000002 => |Sx|<=1.
mace3: |Sz| < (3.1+1e-8+1.5)/2.5 < 1.840000005 => |Sz|<=1;
|Sy| < (3.1+1e-8+1.3+.3)/2.3 < 2.043478266 => |Sy|<=2;
|Sx| < (3.1+1e-8+.9+.7*2+.2)/2.4 < 2.333333338 => |Sx|<=2.
gauge: translation adds <=(3,2,0), so |S|<=(5,4,1).
Direct gauge coordinates also give |Sz|<=1,
|Sy| < (3.1+1e-8+5.4+.3)/2.3 < 3.826086962 => |Sy|<=3,
|Sx| < (3.1+1e-8+6.7+.7*3+.2)/2.4 < 5.041666672 => |Sx|<=5.
deform: |Sz|<=1;
|Sy| < (3.1+1e-8+1.3+.3)/2.35 < 2.000000005 => |Sy|<=2;
|Sx| < (3.1+1e-8+.9+.8*2+.2)/2.4 < 2.416666672 => |Sx|<=2.
displaced: the maximum spans remain (.9,1.3,1.5), so the mace3 bounds apply.
```

For cubic empty/self/multi4, no cross terms occur. Self at R=2.1 gives
(2.1+1e-8)/2 < 1.050000006 on all axes, hence |Sk|<=1. The strict norm
admits exactly the six axial images; R=1.9 gives |Sk|=0 and thus no self
images. Disabling periodic self removes all these tuples. Multi4 has spans
(.6,0,0): x < (1.5+1e-8+.6)/2 < 1.050000006 and y,z < .750000006,
hence bounds (1,0,0) and exactly the four declared directed tuples. Empty
input has no pairs; the same cell/cutoff is valid. Boundary has maximum
x span next_up(2), so x < (2+1e-8+next_up(2))/8 < .500000002 and
y,z < .250000002; all shifts are zero. Its equality/outside positions are
excluded and its next_down position is included by the specified f64 norm.
Removing PBC axes only reduces these sets, so all eight masks inherit the
same certificates. Every component is strictly below 16.

For the later displaced finite-difference suite, moving any Cartesian
coordinate by +/-1e-5 increases any span by at most 2e-5. Substituting spans
(.90002,1.30002,1.50002) gives z<1.840008005, y<2.043486962,
x<2.333341672, retaining (2,2,1). This certificate does not certify a model
or stand in for executing the later model gate.

The triangular argument also controls rounding: start with z, which has one
nonzero diagonal product and no cancellation, to bound even an arbitrary
i64-to-f64 cast; then y and x. A passing finite f64 norm therefore forces
the small integer ranges above, where integer conversion is exact. With
these coordinates, basis entries <=49 and resulting shifts <=5, each
component has at most seven subtraction/product/add operations with absolute
intermediate sum below 1000. A conservative 32*EPS*1000 < 7.2e-12 bound,
including norm classification and the tiny norm-to-radius allowance, is
well below 1e-8. No subnormal fixture coordinates/cutoffs are used. The
search-only cushion covers arithmetic error and boundary-margin exclusion;
no epsilon is used in the oracle membership predicate. The executable
triangular certificate test independently substitutes every fixture's actual
constants and checks the documented integer limits against oracle output.

## Official gate and limitations

`periodic_graph_official::official_ase_images` first checks every official
fixture's cube cutoff margin >=1e-8 and exact Rust/cube tuple equality,
then calls the existing explicitly selected interpreter with actual Rust
edges. ASE 3.29.0 / NumPy 2.5.4 `neighbor_list('ijSD', ...)` supplies the
reference tuples and independent displacement vectors. The helper fails on
pin drift, missing packages, changed tuples, or vector disagreement at
atol=rtol=2e-12. The driver requires response fields and reports bounded
stderr on failures. Each response records package/Python/platform identity,
neighborlist.py path/SHA256, request and canonical fixture SHA256, fixture-file
SHA256 and helper SHA256. An ignored test merely compiling is not acceptance.
These finite fixtures establish G1 coverage only, not model/MD parity or
performance. The external gate must be explicitly executed and retained.
