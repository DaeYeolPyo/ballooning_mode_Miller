import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg, sparse
from scipy.interpolate import PchipInterpolator
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse.linalg import eigsh


class Sturm_Liouville:
    """Solve (g X')' + c X = lambda f X using P2 finite elements."""

    def __init__(
        self,
        thetaGrid,
        g,
        c,
        f,
        bcType="dirichlet",
        numEigenvalues=None,
    ):
        arrays = []

        for name, value in zip(
            ("thetaGrid", "g", "c", "f"),
            (thetaGrid, g, c, f),
        ):
            arr = np.asarray(value, dtype=float)

            if arr.ndim not in (1, 2) or (
                arr.ndim == 2 and 1 not in arr.shape
            ):
                raise ValueError(f"{name} must be a vector.")

            arr = arr.reshape(-1)

            if not np.all(np.isfinite(arr)):
                raise ValueError(f"{name} must contain finite values.")

            arrays.append(arr.copy())

        self.thetaGrid, self.g, self.c, self.f = arrays
        self.ngrid = self.thetaGrid.size

        if self.ngrid < 4:
            raise ValueError("At least four grid points are required.")

        if any(arr.size != self.ngrid for arr in arrays[1:]):
            raise ValueError(
                "thetaGrid, g, c, and f must have the same length."
            )

        if np.any(np.diff(self.thetaGrid) <= 0):
            raise ValueError("thetaGrid must be strictly increasing.")

        self.bcType = (
            "dirichlet"
            if bcType is None or bcType == ""
            else bcType
        )

        self.construct_matrix()
        self.calculate_eigenvalue(numEigenvalues)

    @staticmethod
    def GauLeg_points(nq):
        """Gauss-Legendre quadrature points and weights on [-1, 1]."""
        if not isinstance(nq, (int, np.integer)) or nq < 1:
            raise ValueError("nq must be a positive integer.")

        return np.polynomial.legendre.leggauss(nq)

    def construct_matrix(self):
        nelm = self.ngrid - 1
        self.ndof = 2 * nelm + 1

        # Global P2 DOFs: endpoints and element midpoints.
        self.thetaDof = np.empty(self.ndof)
        self.thetaDof[::2] = self.thetaGrid
        self.thetaDof[1::2] = 0.5 * (
            self.thetaGrid[:-1] + self.thetaGrid[1:]
        )

        xi, w = self.GauLeg_points(4)

        # Reference P2 basis functions and derivatives.
        N = np.column_stack(
            (
                0.5 * xi * (xi - 1.0),
                1.0 - xi**2,
                0.5 * xi * (xi + 1.0),
            )
        )
        dN = np.column_stack(
            (
                xi - 0.5,
                -2.0 * xi,
                xi + 0.5,
            )
        )

        # Element Jacobians and physical quadrature points.
        J = 0.5 * np.diff(self.thetaGrid)
        midpoint = self.thetaDof[1::2]
        theta_q = midpoint[:, None] + J[:, None] * xi

        gq, cq, fq = (
            PchipInterpolator(self.thetaGrid, values)(theta_q)
            for values in (self.g, self.c, self.f)
        )

        # Element matrices.
        # G_ij = integral g * N_i' * N_j' dtheta
        # C_ij = integral c * N_i  * N_j  dtheta
        # F_ij = integral f * N_i  * N_j  dtheta
        Ge = np.einsum(
            "eq,qi,qj->eij",
            gq * w / J[:, None],
            dN,
            dN,
        )
        Ce = np.einsum(
            "eq,qi,qj->eij",
            cq * w * J[:, None],
            N,
            N,
        )
        Fe = np.einsum(
            "eq,qi,qj->eij",
            fq * w * J[:, None],
            N,
            N,
        )

        connectivity = (
            2 * np.arange(nelm)[:, None] + np.arange(3)
        )

        rows = np.broadcast_to(
            connectivity[:, :, None], (nelm, 3, 3)
        ).ravel()
        cols = np.broadcast_to(
            connectivity[:, None, :], (nelm, 3, 3)
        ).ravel()

        shape = (self.ndof, self.ndof)

        # COO -> CSR sums contributions at shared global DOFs.
        self.Gmat, self.Cmat, self.Fmat = (
            sparse.coo_matrix(
                (local.ravel(), (rows, cols)),
                shape=shape,
            ).tocsr()
            for local in (Ge, Ce, Fe)
        )

        return (
            self.Gmat,
            self.Cmat,
            self.Fmat,
            self.thetaDof,
        )

    @staticmethod
    def _bracket_largest_eigenvalue(A, B):
        """
        Find an upper bound for the largest generalized eigenvalue.

        For SPD B, sigma*B - A is SPD exactly when sigma exceeds
        the largest generalized eigenvalue.
        """
        # Reorder to keep the bandwidth small, including periodic BCs.
        permutation = reverse_cuthill_mckee(
            (abs(A) + abs(B)).tocsr(),
            symmetric_mode=True,
        )

        matrices = [
            matrix[permutation, :][:, permutation].tocoo()
            for matrix in (A, B)
        ]

        bandwidth = max(
            int(
                np.max(
                    np.abs(matrix.row - matrix.col),
                    initial=0,
                )
            )
            for matrix in matrices
        )

        bands = []

        for matrix in matrices:
            band = np.zeros((bandwidth + 1, A.shape[0]))
            lower = matrix.row >= matrix.col

            band[
                matrix.row[lower] - matrix.col[lower],
                matrix.col[lower],
            ] = matrix.data[lower]

            bands.append(band)

        Aband, Bband = bands

        def positive_definite(band):
            try:
                linalg.cholesky_banded(
                    band,
                    lower=True,
                    check_finite=False,
                )
                return True
            except linalg.LinAlgError:
                return False

        if not positive_definite(Bband):
            raise ValueError(
                "The mass matrix must be positive definite."
            )

        step = 1.0e-6

        if positive_definite(-Aband):
            upper, lower = 0.0, -step

            for _ in range(80):
                if not positive_definite(
                    lower * Bband - Aband
                ):
                    break

                lower *= 2.0
            else:
                raise RuntimeError(
                    "Could not bracket the largest eigenvalue."
                )

        else:
            lower, upper = 0.0, step

            for _ in range(80):
                if positive_definite(
                    upper * Bband - Aband
                ):
                    break

                upper *= 2.0
            else:
                raise RuntimeError(
                    "Could not bracket the largest eigenvalue."
                )

        for _ in range(80):
            if upper - lower <= (
                1.0e-10 * max(1.0, abs(upper))
            ):
                break

            middle = 0.5 * (lower + upper)

            if positive_definite(
                middle * Bband - Aband
            ):
                upper = middle
            else:
                lower = middle

        return upper

    def calculate_eigenvalue(self, numEigenvalues=None):
        """
        Solve (C - G) X = lambda F X.

        numEigenvalues=None or np.inf: all eigenvalues.
        Positive integer: that many largest algebraic eigenvalues.
        """
        if numEigenvalues is None:
            requested = np.inf
        else:
            requested = float(numEigenvalues)

            if (
                np.isnan(requested)
                or requested < 1
                or (
                    np.isfinite(requested)
                    and requested != np.floor(requested)
                )
            ):
                raise ValueError(
                    "numEigenvalues must be a positive integer or Inf."
                )

        self.Amat = self.Cmat - self.Gmat
        self.Bmat = self.Fmat

        ndof = self.ndof
        bc = self.bcType.lower()

        # X_full = P @ X_reduced
        if bc in ("natural", "neumann"):
            # g X' = 0 at both endpoints.
            P = sparse.eye(ndof, format="csr")

        elif bc == "dirichlet":
            # X(left) = X(right) = 0.
            P = sparse.eye(ndof, format="csr")[:, 1:-1]

        elif bc == "natural-dirichlet":
            # Natural at left; X(right) = 0.
            P = sparse.eye(ndof, format="csr")[:, :-1]

        elif bc == "periodic":
            # Identify the left and right endpoint DOFs.
            rows = np.arange(ndof)
            cols = np.r_[np.arange(ndof - 1), 0]

            P = sparse.coo_matrix(
                (np.ones(ndof), (rows, cols)),
                shape=(ndof, ndof - 1),
            ).tocsr()

        else:
            raise ValueError(
                f"Unknown boundary condition: {self.bcType}"
            )

        A = (P.T @ self.Amat @ P).tocsr()
        B = (P.T @ self.Bmat @ P).tocsr()

        size = A.shape[0]
        scale = np.ones(size)

        if np.isfinite(requested) and requested < size - 1:
            for matrix in (A, B):
                symmetry_error = sparse.linalg.norm(
                    matrix - matrix.T
                )
                tolerance = 1.0e-12 * max(
                    sparse.linalg.norm(matrix),
                    np.finfo(float).tiny,
                )

                if symmetry_error > tolerance:
                    raise ValueError(
                        "Sparse eigenvalue solving requires "
                        "symmetric matrices."
                    )

            if np.any(B.diagonal() <= 0):
                raise ValueError(
                    "The mass matrix must have a positive diagonal."
                )

            # Diagonal scaling, as in the MATLAB solver.
            scale = 1.0 / np.sqrt(B.diagonal())
            D = sparse.diags(scale)

            A = (D @ A @ D).tocsr()
            B = (D @ B @ D).tocsr()

            options = dict(
                M=B,
                k=int(requested),
                tol=1.0e-10,
                maxiter=4000,
                v0=np.ones(size) / np.sqrt(size),
            )

            if requested == 1:
                upper = self._bracket_largest_eigenvalue(A, B)
                shift = upper + 1.0e-6 * max(
                    1.0, abs(upper)
                )

                values, vectors = eigsh(
                    A,
                    sigma=shift,
                    which="LM",
                    **options,
                )
            else:
                values, vectors = eigsh(
                    A,
                    which="LA",
                    **options,
                )

        else:
            # Generalized dense eigensolver, matching MATLAB eig(A, B).
            values, vectors = linalg.eig(
                A.toarray(),
                B.toarray(),
            )

        # Preserve the MATLAB convention: descending real eigenvalues.
        order = np.argsort(values.real)[::-1]

        if np.isfinite(requested):
            order = order[:int(requested)]

        self.eigval = values.real[order]
        self.eigfunc = np.real_if_close(
            P @ (scale[:, None] * vectors[:, order])
        )

        return self.eigval, self.eigfunc


def _plot_test(
    eq,
    mode,
    lambdaExact,
    analytic,
    numplot,
    title,
    multiplier=1.0,
):
    mode = np.asarray(mode)

    if mode.size > eq.eigval.size:
        raise ValueError(
            "nMode exceeds the number of available eigenvalues."
        )

    if not 0 <= numplot <= mode.size:
        raise ValueError(
            "numplot must be between 0 and nMode."
        )

    fig, ax = plt.subplots(num=f"{title}: eigenvalues")
    ax.plot(
        mode,
        eq.eigval[:mode.size],
        "bo--",
        label="Numerical",
    )
    ax.plot(
        mode,
        lambdaExact,
        "r-",
        label="Analytic",
    )
    ax.set(
        xlabel="Mode",
        ylabel=r"$\lambda$",
        title=title,
    )
    ax.legend()
    ax.grid(True)

    if numplot:
        fig, axes = plt.subplots(
            numplot,
            1,
            squeeze=False,
            figsize=(8, 2.5 * numplot),
            num=f"{title}: eigenfunctions",
        )

        for i, ax in enumerate(axes[:, 0]):
            Xnum = multiplier * eq.eigfunc[:, i]
            Xana = analytic(i)

            # Normalize using the same Fmat metric as MATLAB.
            Xnum = Xnum / np.sqrt(
                np.real(Xnum.conj() @ (eq.Fmat @ Xnum))
            )
            Xana = Xana / np.sqrt(
                np.real(Xana.conj() @ (eq.Fmat @ Xana))
            )

            # Eigenvectors have an arbitrary overall sign.
            if np.real(
                Xana.conj() @ (eq.Fmat @ Xnum)
            ) < 0:
                Xnum = -Xnum

            ax.plot(
                eq.thetaDof,
                Xnum,
                "b--",
                label="Numerical",
            )
            ax.plot(
                eq.thetaDof,
                Xana,
                "r-",
                label="Analytic",
            )
            ax.set(
                xlabel=r"$\theta$",
                ylabel="X",
                title=f"Mode {mode[i]}",
            )
            ax.legend()
            ax.grid(True)

        fig.tight_layout()

    return eq


def test_assoc_Legendre(
    ngrid=200,
    mOrder=1,
    nMode=100,
    numplot=5,
):
    from scipy.special import lpmv

    if (
        not isinstance(mOrder, (int, np.integer))
        or mOrder < 0
    ):
        raise ValueError(
            "mOrder must be a nonnegative integer."
        )

    thetaGrid = np.linspace(0.0, np.pi, ngrid)
    weight = np.sin(thetaGrid)**(2 * mOrder + 1)

    eq = Sturm_Liouville(
        thetaGrid,
        weight,
        -mOrder * (mOrder + 1) * weight,
        weight,
        bcType="natural",
    )

    mode = mOrder + np.arange(nMode)

    # Physical associated Legendre function:
    # P_l^m(cos(theta)) = sin(theta)^m * X(theta).
    return _plot_test(
        eq,
        mode,
        -mode * (mode + 1),
        lambda i: lpmv(
            mOrder,
            int(mode[i]),
            np.cos(eq.thetaDof),
        ),
        numplot,
        "Associated Legendre",
        multiplier=np.sin(eq.thetaDof)**mOrder,
    )


def test_Legendre(
    ngrid=200,
    nMode=100,
    numplot=5,
):
    from scipy.special import eval_legendre

    thetaGrid = np.linspace(0.0, np.pi, ngrid)
    weight = np.sin(thetaGrid)

    eq = Sturm_Liouville(
        thetaGrid,
        weight,
        np.zeros(ngrid),
        weight,
        bcType="natural",
    )

    mode = np.arange(nMode)

    return _plot_test(
        eq,
        mode,
        -mode * (mode + 1),
        lambda i: eval_legendre(
            int(mode[i]),
            np.cos(eq.thetaDof),
        ),
        numplot,
        "Legendre",
    )


def test_Fourier(
    ngrid=200,
    c0=3.0,
    nMode=30,
    numplot=5,
):
    thetaGrid = np.linspace(
        0.0, 2.0 * np.pi, ngrid
    )

    eq = Sturm_Liouville(
        thetaGrid,
        np.ones(ngrid),
        np.full(ngrid, c0),
        np.ones(ngrid),
        bcType="periodic",
    )

    mode = np.arange(1, nMode + 1)
    harmonics = (np.arange(nMode) + 1) // 2

    def analytic(i):
        m = harmonics[i]

        if m == 0:
            return np.ones(eq.ndof)

        # Degenerate Fourier modes may be arbitrary combinations
        # of cos(m*theta) and sin(m*theta).
        Phi = np.column_stack(
            (
                np.cos(m * eq.thetaDof),
                np.sin(m * eq.thetaDof),
            )
        )

        coefficients = np.linalg.solve(
            Phi.T @ (eq.Fmat @ Phi),
            Phi.T @ (
                eq.Fmat @ eq.eigfunc[:, i]
            ),
        )

        return Phi @ coefficients

    return _plot_test(
        eq,
        mode,
        c0 - harmonics**2,
        analytic,
        numplot,
        "Fourier",
    )


def test_Bessel(
    ngrid=200,
    mOrder=1,
    R=2.0 * np.pi,
    nMode=100,
    numplot=5,
):
    from scipy.special import jn_zeros, jv

    if (
        not isinstance(mOrder, (int, np.integer))
        or mOrder < 0
    ):
        raise ValueError(
            "mOrder must be a nonnegative integer."
        )

    if not np.isfinite(R) or R <= 0:
        raise ValueError(
            "R must be positive and finite."
        )

    thetaGrid = np.linspace(0.0, R, ngrid)
    weight = thetaGrid**(2 * mOrder + 1)

    eq = Sturm_Liouville(
        thetaGrid,
        weight,
        np.zeros(ngrid),
        weight,
        bcType="natural-dirichlet",
    )

    mode = np.arange(1, nMode + 1)
    zeros = jn_zeros(mOrder, nMode)

    # Physical Bessel function:
    # J_m(j_mn*theta/R) = theta^m * X(theta).
    return _plot_test(
        eq,
        mode,
        -(zeros / R)**2,
        lambda i: jv(
            mOrder,
            zeros[i] * eq.thetaDof / R,
        ),
        numplot,
        "Bessel",
        multiplier=eq.thetaDof**mOrder,
    )


def test_inhom(
    ngrid=200,
    c0=0.0,
    nMode=100,
    numplot=5,
):
    thetaGrid = np.linspace(0.0, 1.0, ngrid)
    weight = 1.0 / (1.0 + thetaGrid)

    eq = Sturm_Liouville(
        thetaGrid,
        1.0 + thetaGrid,
        c0 * weight,
        weight,
        bcType="dirichlet",
    )

    mode = np.arange(1, nMode + 1)

    return _plot_test(
        eq,
        mode,
        c0 - (mode * np.pi / np.log(2.0))**2,
        lambda i: np.sin(
            mode[i]
            * np.pi
            * np.log1p(eq.thetaDof)
            / np.log(2.0)
        ),
        numplot,
        "Inhomogeneous",
    )


if __name__ == "__main__":
    eq = test_Bessel(
        ngrid=200,
        mOrder=1,
        nMode=100,
        numplot=5,
    )

    plt.show()
