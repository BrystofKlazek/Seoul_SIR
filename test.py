import time
import numpy as np
import SIR_diffusion_wrap as sir
import matplotlib.pyplot as plt

import numpy as np

def grid_laplacian(nx, ny, diagonals=False, w_card=1, w_diag=1):
    n = nx * ny
    L = np.zeros((n, n), dtype=np.float64)

    def idx(i, j):
        return i * ny + j

    for i in range(nx):
        for j in range(ny):
            p = idx(i, j)
            deg = 0.0

            if 0 <= i-1 < nx:
                q = idx(i-1, j)
                L[p, q] -= w_card
                deg += w_card

            if 0 <= i+1 < nx:
                q = idx(i+1, j)
                L[p, q] -= w_card
                deg += w_card

            if 0 <= j-1 < ny:
                q = idx(i, j-1)
                L[p, q] -= w_card
                deg += w_card

            if 0 <= j+1 < ny:
                q = idx(i, j+1)
                L[p, q] -= w_card
                deg += w_card

            if diagonals:
                if 0 <= i-1 < nx and 0 <= j-1 < ny:
                    q = idx(i-1, j-1)
                    L[p, q] -= w_diag
                    deg += w_diag

                if 0 <= i-1 < nx and 0 <= j+1 < ny:
                    q = idx(i-1, j+1)
                    L[p, q] -= w_diag
                    deg += w_diag

                if 0 <= i+1 < nx and 0 <= j-1 < ny:
                    q = idx(i+1, j-1)
                    L[p, q] -= w_diag
                    deg += w_diag

                if 0 <= i+1 < nx and 0 <= j+1 < ny:
                    q = idx(i+1, j+1)
                    L[p, q] -= w_diag
                    deg += w_diag

            L[p, p] = deg

    return L

def main():
    nx = ny = 50
    dx = dy = 1.0
    dt = 0.01
    steps = 200
    snapshots = 0  # set to 0; wrapper ignores it internally for now, I will add snapshots later

    print("Building test data")
    x0 = np.zeros((nx, ny), dtype=np.float64)
    x0[nx//2, ny//2] = 1.0

    out_grid = sir.euler_solve_normal(x0, dx, dy, dt, steps, snapshots)
    print(f"min={out_grid.min():.3e}, max={out_grid.max():.3e}")

    print("Building graph Laplacian ")
    L = grid_laplacian(nx, ny)  
    out_graph = sir.euler_solve_graph(L, x0, dt, steps, snapshots)
    print(f"min={out_graph.min():.3e}, max={out_graph.max():.3e}")

    out_graph_grid = out_graph.reshape(nx, ny)

    vmin = min(out_grid.min(), out_graph_grid.min())
    vmax = max(out_grid.max(), out_graph_grid.max())

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), constrained_layout=True)

    im0 = axes[0].imshow(out_grid, origin="lower", cmap="viridis",
                         vmin=vmin, vmax=vmax, interpolation="nearest")
    axes[0].set_title("Grid solver (Euler)")
    axes[0].set_xlabel("j"); axes[0].set_ylabel("i")

    im1 = axes[1].imshow(out_graph_grid, origin="lower", cmap="viridis",
                         vmin=vmin, vmax=vmax, interpolation="nearest")
    axes[1].set_title("Graph solver reshaped")
    axes[1].set_xlabel("j"); axes[1].set_ylabel("i")

    diff = out_graph_grid - out_grid
    im2 = axes[2].imshow(diff, origin="lower", cmap="coolwarm",
                         interpolation="nearest")
    axes[2].set_title("Difference (graph - grid)")
    axes[2].set_xlabel("j"); axes[2].set_ylabel("i")

    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04).set_label("u")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04).set_label("u")
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04).set_label("Δu")

    plt.show()

if __name__ == "__main__":
    main()

