# Geometry stage

This folder defines and meshes only the prescribed last closed flux surface
(LCFS). It does not prescribe any interior flux surface.

## Files

- `make_boundary.m` supports the original harmonic C-shape, a localized
  outboard-midplane indentation, a two-arc pointed boomerang, and a
  third-harmonic cosine boundary, a Miller D-shape, and ordered tabulated
  LCFS points such as a GEQDSK boundary. It rejects non-positive-radius and
  self-intersecting boundaries and reports tip, endpoint, bounding-box
  elongation, and other geometric diagnostics.
- `generate_mesh.m` fills a convex or concave LCFS with a constrained,
  counterclockwise P1 triangular mesh. It returns ordered boundary nodes for
  the later fixed Dirichlet condition.
- `verify_geometry.m` checks the analytic and tabulated parameterizations.

From the `GS_fixed` directory, run:

```matlab
addpath('geometry')
report = verify_geometry(true);
```

To create a localized extreme C-shaped case whose tip lies inside `A=3`:

```matlab
params = struct('type','localized-c', ...
                'A',3.0,'B',0.9,'D',1.1,'beta',0.9, ...
                'G',1.3298,'H',0.0,'nBoundary',400);
geom = make_boundary(params);
mesh = generate_mesh(geom, struct('targetH',0.05));
```

The localized formula is

```text
R = A + B sin(theta) - D exp(beta (sin(theta)-1)),
Z = G cos(theta) - H sin(2 theta).
```

At the outboard midplane, `R_tip=A+B-D`; therefore `D>B` places the
indentation tip inside the nominal major radius `A`. `beta` controls angular
localization: larger values make a narrower notch while leaving the inboard
boundary nearly unchanged. Very narrow, deep notches can force multiple
O/X-points even when the polygon is valid. For `D=1.1` and elongation about
3, `beta=0.9` is the verified concave, single-axis reference; increase beta
only with mesh and topology scans.

`targetH` controls the nominal interior-node spacing. Boundary resolution is
controlled independently by `nBoundary`; this separation is intentional.

For an imported ordered LCFS, use:

```matlab
params = struct('type','tabulated','R',rbbbs,'Z',zbbbs, ...
                'nBoundary',numel(rbbbs));
geom = make_boundary(params);
```

When `nBoundary` equals the number of distinct input points, those points are
preserved exactly. Otherwise the closed curve is resampled uniformly in
arclength.

## Pointed boomerang C-shape

For an inward midplane indentation and outward-swept pointed ends, use:

```matlab
params = struct('type','boomerang-c','A',3.0, ...
                'Rback',2.0,'Rtip',2.65, ...
                'Rend',3.8,'Zend',2.7, ...
                'pInner',1.6,'pOuter',1.4, ...
                'nBoundary',400);
geom = make_boundary(params);
```

With `u=|Z|/Zend`, the two arcs are

```text
R_inner = Rtip  + (Rend-Rtip ) u^pInner,
R_back  = Rback + (Rend-Rback) u^pOuter.
```

They meet at `(Rend,+/-Zend)`, so the upper and lower ends are genuine
polygon corners rather than rounded harmonic extrema. At the midplane the
plasma spans `Rback <= R <= Rtip`. The reference values give `Rtip=2.65<A=3`,
endpoints at `Rend=3.8`, and bounding-box elongation
`2*Zend/(Rend-Rback)=3`. `pInner` and `pOuter` control the sweep of the two
arcs; invalid combinations that make them cross are rejected. Exact pointed
ends concentrate mesh gradients, so increase both `nBoundary` and the
interior resolution (`targetH` smaller) for quantitative edge studies.

## Third-harmonic cosine C-shape

The proposed Fourier-like boundary is available as `type='cosine-c'`:

```matlab
params = struct('type','cosine-c', ...
                'R0',3.0,'A',0.8,'B',-1.0,'C',0.0, ...
                'kappa',4.0,'nBoundary',400);
geom = make_boundary(params);
```

```text
R = R0 + A cos(theta) + B cos(2 theta) + C cos(3 theta),
Z = kappa A sin(theta).
```

The exact vertical height is `2*kappa*A`. The two midplane radii are
`R0+B+(A+C)` and `R0+B-(A+C)`, while the upper/lower point is at
`(R0-B,+/-kappa*A)`. Thus negative `B` moves both midplane crossings inward
and the upper/lower points outward. The reported `nominalElongation` is
`kappa`; `boundingBoxElongation` is computed from the actual radial extrema
and generally differs from `kappa` once `B` or `C` is nonzero.

The upper/lower points become exact cusps when `C=A/3`. Such a zero-opening
tip produces very low-quality triangles and is not recommended for a P1 FEM
solve without local tip regularization. Values near, but not exactly at,
`A/3` give a sharp rounded end.

`geom.outboardMidplaneCurvatureRadius` reports the local poloidal-plane
curvature radius at the outboard midplane. For the usual `A+C>0` branch it is

```text
R_curv,tip = (kappa A)^2 / |A + 4 B + 9 C|.
```
