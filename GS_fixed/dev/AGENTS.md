# AGENTS.md

## Project overview
This repository develops a MATLAB-based high-order finite-element
Grad-Shafranov solver for fixed-boundary tokamak equilibria.

The current main discretization is P3 Lagrange triangular FEM.

## Before modifying code
Read:
- docs/CODEX_CONTEXT.md
- docs/NUMERICAL_DESIGN.md

Do not redesign the numerical formulation without checking these files.

## MATLAB conventions
- Use MATLAB functions rather than scripts for reusable components.
- Prefer vectorized calculations when readability is maintained.
- Use sparse matrices for global FEM matrices.
- Preserve existing function signatures unless necessary.
- Add input validation where appropriate.

## FEM conventions
Refernce triangle:

lambda1 = 1 - xi - eta
lambda2 = xi
lambda3 = eta

P3 elements contain 10 local nodes.

Use the node ordering defined in:
docs/CODEX_CONTEXT.md

Do not silently change local node ordering.

## Quadrature
Use Gauss_Dunavant_quadrature.m for triangle integration.

Choose a quadrature degree sufficient for the polynomial degree
of the integrand.

## Verification
Whenever changing FEM routines:
1. Check partition of unity.
2. Check Kronecker-delta property of basis functions.
3. Check symmetry of the stiffness matrix.
4. Run existing verification scripts.
5. Avoid modifying tests merely to make them pass.

## Coding philosophy
Prioritize numerical correctness over premature optimization.

When vectorizing code, keep the mathematical correspondence
with the FEM formulation understandable.

## Current project context
See:
docs/CODEX_CONTEXT.md