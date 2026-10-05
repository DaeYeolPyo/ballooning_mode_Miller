# Computed original C-shape dense map

Updated 2026-10-02 21:55:10.

Attempted **374/374** requested points; 374 actual beta values; 374 complete stability scans.

Stable 147; unstable 227; marginal 0; unknown 0.

Actual volume beta range: 0.065865% to 13.258216%.

Maximum accepted GS relative residual: 4.76892e-07; maximum accepted q relative error: 0.00199223.

The map uses actual beta = 2 mu0 integral(p dV)/integral(B^2 dV), including both poloidal and toroidal field.

The stability classification samples 11 surfaces and 8 phases per equilibrium.

Numerical resolution varies only where recovery is needed. Each CSV row records the actual GS/q and surface/Fourier controls; the MAT point records full solverConfiguration and ballooning.options. Acceptance tolerances and physical profiles are unchanged.

Refined 17/17 transition brackets in 17 columns, using 85 extra equilibrium/stability evaluations.

| I | d_out/(2 a_g) | transition | beta_av bracket [%] | converged |
|---:|---:|:---|:---|:---|
| 0.600 | 0.078988 | stable to unstable | 3.372583 to 3.385122 | 1 |
| 0.625 | 0.088402 | stable to unstable | 3.323556 to 3.335871 | 1 |
| 0.650 | 0.098201 | stable to unstable | 3.261673 to 3.273986 | 1 |
| 0.675 | 0.108372 | stable to unstable | 3.211797 to 3.224100 | 1 |
| 0.700 | 0.118901 | stable to unstable | 3.160987 to 3.173226 | 1 |
| 0.725 | 0.129776 | stable to unstable | 3.122950 to 3.135056 | 1 |
| 0.750 | 0.140987 | stable to unstable | 3.073360 to 3.085474 | 1 |
| 0.775 | 0.152523 | stable to unstable | 3.022749 to 3.034729 | 1 |
| 0.800 | 0.164373 | stable to unstable | 2.972345 to 2.984261 | 1 |
| 0.825 | 0.176530 | stable to unstable | 2.920485 to 2.932313 | 1 |
| 0.850 | 0.188985 | stable to unstable | 2.870731 to 2.882516 | 1 |
| 0.875 | 0.201730 | stable to unstable | 2.807213 to 2.818914 | 1 |
| 0.900 | 0.214757 | stable to unstable | 2.758997 to 2.770599 | 1 |
| 0.925 | 0.228060 | stable to unstable | 2.696889 to 2.708371 | 1 |
| 0.950 | 0.241631 | stable to unstable | 2.634640 to 2.646082 | 1 |
| 0.975 | 0.255466 | stable to unstable | 2.541520 to 2.552899 | 1 |
| 1.000 | 0.269556 | stable to unstable | 2.467453 to 2.478669 | 1 |

These are the current solver results, not an independent validation of the high-beta ballooning formulation. See ../../indented_cshape_dense_README.md for equations, resolution, and resume commands.
