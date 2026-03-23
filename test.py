import time
import numpy as np
import SIR_diffusion_wrap as sir
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# ----- CSR Laplacian for an nx x ny grid -----
def grid_laplacian_csr(nx, ny, diagonals=False, w_card=1.0, w_diag=1.0, mirror=True):
    n = nx * ny
    indptr = [0]
    indices = []
    data = []

    def idx(i, j): return i * ny + j

    for i in range(nx):
        for j in range(ny):
            p = idx(i, j)
            deg = 0.0

            def add_edge(q, w):
                nonlocal deg
                indices.append(q)
                data.append(-w)
                deg += w

            has_im = (i-1 >= 0)
            has_ip = (i+1 < nx)
            if has_im: add_edge(idx(i-1, j), w_card)
            if has_ip: add_edge(idx(i+1, j), w_card)
            if mirror:
                if not has_im and has_ip:  
                    add_edge(idx(i+1, j), w_card)  
                if not has_ip and has_im:  
                    add_edge(idx(i-1, j), w_card)  

            # y-direction
            has_jm = (j-1 >= 0)
            has_jp = (j+1 < ny)
            if has_jm: add_edge(idx(i, j-1), w_card)
            if has_jp: add_edge(idx(i, j+1), w_card)
            if mirror:
                if not has_jm and has_jp:  
                    add_edge(idx(i, j+1), w_card)
                if not has_jp and has_jm:  
                    add_edge(idx(i, j-1), w_card)

            if diagonals:
                for di, dj in [(-1,-1), (-1,1), (1,-1), (1,1)]:
                    ii, jj = i+di, j+dj
                    if 0 <= ii < nx and 0 <= jj < ny:
                        add_edge(idx(ii, jj), w_diag)

            indices.append(p)
            data.append(deg)

            indptr.append(len(indices))

    return (np.asarray(indptr, dtype=np.int32),
            np.asarray(indices, dtype=np.int32),
            np.asarray(data, dtype=np.float64))

def sir_rhs(state, out, beta=0.8, gamma=0.2):
    S, I, R = state[0], state[1], state[2]
    dS, dI, dR = out[0], out[1], out[2]
    inf = beta * S * I
    dS[...] = -inf
    dI[...] =  inf - gamma * I
    dR[...] =  gamma * I

def main():
    nx = ny = 50
    dx = dy = 1.0
    dt = 0.02
    steps = 1000
    snapshots = 50
    beta, gamma = 0.8, 0.2

    S0 = np.ones((nx, ny), np.float64)
    I0 = np.zeros((nx, ny), np.float64); I0[nx//2, ny//2] = 1e-2
    R0 = np.zeros((nx, ny), np.float64)
    x0_grid = np.stack([S0, I0, R0], axis=0)

    print("Solving on grid (Euler, reaction-diffusion SIR)...")
    t0 = time.time()
    out_grid = sir.euler_solve_normal_rd(
        x0_grid, dx, dy, dt, steps, snapshots,
        rhs_fn=sir_rhs, params={"beta": beta, "gamma": gamma}
    )
    print(f"grid done in {time.time()-t0:.2f}s, snaps={out_grid.shape[0]}")  # (T, 3, nx, ny)

    print("Building CSR Laplacian…")
    indptr, indices, data = grid_laplacian_csr(nx, ny, diagonals=False, mirror=True)

    # The CSR solver expects a list of arrays per time step (time-varying graph).
    # For a static graph, repeat the same arrays for every step.
    indptr_list  = [indptr]  * steps
    indices_list = [indices] * steps
    data_list    = [data]    * steps

    x0_graph = np.stack([S0.ravel(), I0.ravel(), R0.ravel()], axis=0)

    print("Solving on graph (CSR, Euler, reaction-diffusion SIR)…")
    t0 = time.time()
    out_graph = sir.euler_solve_graph_csr_rd(
        indptr_list, indices_list, data_list, x0_graph,
        dt, steps, snapshots,
        rhs_fn=sir_rhs, params={"beta": beta, "gamma": gamma}
    )
    print(f"graph done in {time.time()-t0:.2f}s, snaps={out_graph.shape[0]}")  # (T, 3, n)

    # Extract S, I, R
    Sg, Ig, Rg = out_grid[:,0], out_grid[:,1], out_grid[:,2]              # (T, nx, ny)
    Sx = out_graph[:,0,:].reshape(-1, nx, ny)
    Ix = out_graph[:,1,:].reshape(-1, nx, ny)
    Rx = out_graph[:,2,:].reshape(-1, nx, ny)

    T = Sg.shape[0]

    # Color ranges per compartment (helps readability)
    vS = (Sg.min(), Sg.max())
    vI = (Ig.min(), Ig.max())
    vR = (Rg.min(), Rg.max())

    fig, axes = plt.subplots(3, 3, figsize=(14, 12), constrained_layout=True)

    # S row
    im_Sg  = axes[0,0].imshow(Sg[0], origin="lower", cmap="viridis", vmin=vS[0], vmax=vS[1], interpolation="nearest")
    im_Sx  = axes[0,1].imshow(Sx[0], origin="lower", cmap="viridis", vmin=vS[0], vmax=vS[1], interpolation="nearest")
    im_Sd  = axes[0,2].imshow(Sx[0]-Sg[0], origin="lower", cmap="coolwarm", interpolation="nearest")
    axes[0,0].set_title("S (grid)"); axes[0,1].set_title("S (graph CSR)"); axes[0,2].set_title("S diff")

    # I row
    im_Ig  = axes[1,0].imshow(Ig[0], origin="lower", cmap="viridis", vmin=vI[0], vmax=vI[1], interpolation="nearest")
    im_Ix  = axes[1,1].imshow(Ix[0], origin="lower", cmap="viridis", vmin=vI[0], vmax=vI[1], interpolation="nearest")
    im_Id  = axes[1,2].imshow(Ix[0]-Ig[0], origin="lower", cmap="coolwarm", interpolation="nearest")
    axes[1,0].set_title("I (grid)"); axes[1,1].set_title("I (graph CSR)"); axes[1,2].set_title("I diff")

    # R row
    im_Rg  = axes[2,0].imshow(Rg[0], origin="lower", cmap="viridis", vmin=vR[0], vmax=vR[1], interpolation="nearest")
    im_Rx  = axes[2,1].imshow(Rx[0], origin="lower", cmap="viridis", vmin=vR[0], vmax=vR[1], interpolation="nearest")
    im_Rd  = axes[2,2].imshow(Rx[0]-Rg[0], origin="lower", cmap="coolwarm", interpolation="nearest")
    axes[2,0].set_title("R (grid)"); axes[2,1].set_title("R (graph CSR)"); axes[2,2].set_title("R diff")

    for i in range(3):
        for j in range(3):
            axes[i,j].set_xlabel("j")
        axes[i,0].set_ylabel("i", rotation=0, labelpad=10)

    fig.colorbar(im_Sg, ax=axes[0,0], fraction=0.046, pad=0.04).ax.set_title("S")
    fig.colorbar(im_Sx, ax=axes[0,1], fraction=0.046, pad=0.04).ax.set_title("S")
    fig.colorbar(im_Sd, ax=axes[0,2], fraction=0.046, pad=0.04).ax.set_title("ΔS")

    fig.colorbar(im_Ig, ax=axes[1,0], fraction=0.046, pad=0.04).ax.set_title("I")
    fig.colorbar(im_Ix, ax=axes[1,1], fraction=0.046, pad=0.04).ax.set_title("I")
    fig.colorbar(im_Id, ax=axes[1,2], fraction=0.046, pad=0.04).ax.set_title("ΔI")

    fig.colorbar(im_Rg, ax=axes[2,0], fraction=0.046, pad=0.04).ax.set_title("R")
    fig.colorbar(im_Rx, ax=axes[2,1], fraction=0.046, pad=0.04).ax.set_title("R")
    fig.colorbar(im_Rd, ax=axes[2,2], fraction=0.046, pad=0.04).ax.set_title("ΔR")

    def update(frame):
        im_Sg.set_data(Sg[frame]); im_Sx.set_data(Sx[frame]); im_Sd.set_data(Sx[frame]-Sg[frame])
        im_Ig.set_data(Ig[frame]); im_Ix.set_data(Ix[frame]); im_Id.set_data(Ix[frame]-Ig[frame])
        im_Rg.set_data(Rg[frame]); im_Rx.set_data(Rx[frame]); im_Rd.set_data(Rx[frame]-Rg[frame])
        return (im_Sg, im_Sx, im_Sd, im_Ig, im_Ix, im_Id, im_Rg, im_Rx, im_Rd)

    anim = FuncAnimation(fig, update, frames=T, interval=100, blit=False)
    PillowWriter(fps=10).setup(fig, "comparison_SIR_all.gif", dpi=100)
    anim.save("comparison_SIR_all.gif", writer=PillowWriter(fps=10))
    print("Saved GIF to comparison_SIR_all.gif")
    plt.show()

if __name__ == "__main__":
    main()
