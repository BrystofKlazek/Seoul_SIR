import networkx as nx
import numpy as np

def build_laplacian_timeseries_for_C(G, hourly_weights,
                                     n_steps, dt, t0=0.0):
    # precompute CSR Laplacians L(t_k) for all time steps

    node_index = nx.get_node_attributes(G, "index")
    nodes = sorted(G.nodes(), key=lambda n: node_index[n])
    n = len(nodes)

    # neighbours for each node: only outgoing edges (successors)
    neighbours_idx = []
    for node in nodes:
        outs = set(G.successors(node))
        neighs = sorted(node_index[nb] for nb in outs)
        neighbours_idx.append(neighs)

    # CSR pattern for the C code
    indptr = [0]
    indices = []
    for i, neighs in enumerate(neighbours_idx):
        indices.extend(neighs)   # off-diagonals 
        indices.append(i)        # diagonal
        indptr.append(len(indices))

    indptr_arr = np.asarray(indptr, dtype=np.int32)
    indices_arr = np.asarray(indices, dtype=np.int32)
    nnz = indices_arr.size

    print(f"Graph has {n} nodes, CSR nnz per step = {nnz}")

    # pattern is reused, only data change
    indptr_list = [indptr_arr] * n_steps
    indices_list = [indices_arr] * n_steps

    data_list = []

    edge_weight_fn = make_edge_weight_fn(hourly_weights, dt)

    for step in range(n_steps):
        data_step = np.empty(nnz, dtype=float)
        pos = 0
        for i, node in enumerate(nodes):
            deg = 0.0
            for j in neighbours_idx[i]:
                nb = nodes[j]

                # directed weight: only node -> nb, no symmetrization
                w = edge_weight_fn(step, node, nb)

                data_step[pos] = -w   # off-diagonal L_ij
                pos += 1
                deg += w              # accumulate outgoing degree

            data_step[pos] = deg      # diagonal L_ii = sum_j w(i -> j)
            pos += 1

        data_list.append(data_step.astype(np.float64, copy=False))

    return indptr_list, indices_list, data_list


