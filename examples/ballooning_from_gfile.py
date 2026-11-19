"""Evaluate one local-geometry ballooning surface over several phases."""
import argparse
import json
import numpy as np

from eq_process import Equilibrium
from salpha import prepare_salpha_equilibrium, evaluate_salpha_point


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gfile")
    parser.add_argument("--psi-n", type=float, default=.5)
    parser.add_argument("--model", choices=["direct", "miller", "cshape", "indented-miller"], default="direct")
    parser.add_argument("--points", type=int, default=501)
    args = parser.parse_args()
    state = prepare_salpha_equilibrium(Equilibrium.from_geqdsk(args.gfile),
        target_psiN=args.psi_n, model=args.model, n_surface_points=512)
    result = evaluate_salpha_point(state, state.shat, state.alpha,
        theta0=np.linspace(-np.pi, np.pi, 8, endpoint=False),
        theta_boundary=8*np.pi, n_balloon_points=args.points)
    if not np.isfinite(result.lambdaMax):
        raise RuntimeError(str(result.failMessage))
    print(json.dumps(dict(psi_n=args.psi_n, lambda_max=float(result.lambdaMax),
        shat=state.shat, alpha=state.alpha,
        interpretation="one sampled surface; positive lambda denotes instability"), indent=2))


if __name__ == "__main__":
    main()
