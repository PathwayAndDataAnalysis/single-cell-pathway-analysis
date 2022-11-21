import numpy
import pandas as pd
import numpy as np
import math
import scipy.stats as stats
from sklearn.utils.extmath import randomized_svd
from scipy.spatial import distance


def get_min_elem(que, dist):
    minVal = math.inf
    minElem = -1
    for i in que:
        if dist[i] < minVal:
            minVal = dist[i]
            minElem = i

    return minElem


def get_shortest_path(mat, sourceInd, targetInd):
    nCells = len(mat)
    dist = [math.inf] * nCells
    prev = [-1] * nCells

    dist[sourceInd] = 0

    que = [*range(0, nCells)]

    while len(que) > 0:
        u = get_min_elem(que, dist)
        if u == targetInd:
            break

        que.remove(u)

        for v in que:
            alt = dist[u] + mat[u][v]
            if alt < dist[v]:
                dist[v] = alt
                prev[v] = u

    return prev


def get_shortest_path_multi_targ(mat, sourceInd, targetInds):
    nCells = len(mat)
    dist = [math.inf] * nCells
    prev = [-1] * nCells

    dist[sourceInd] = 0

    que = [*range(0, nCells)]
    targetInds = set(targetInds)

    while len(targetInds) > 0:
        u = get_min_elem(que, dist)

        if u in targetInds:
            targetInds.remove(u)

        que.remove(u)

        for v in que:
            alt = dist[u] + mat[u][v]
            if alt < dist[v]:
                dist[v] = alt
                prev[v] = u

    return prev


def get_path_as_indices(prev, source, target):
    path = [target]
    while path[0] != source:
        path.insert(0, prev[path[0]])

    return path


def get_path_as_cell_names(path_as_ind, cells):
    path = []
    for i in path_as_ind:
        path.append(cells[i])
    return path


def get_pseudotime(path_as_ind, dist_mat):
    time = [0]
    prevT = 0

    for i in range(1, len(path_as_ind)):
        c1 = path_as_ind[i - 1]
        c2 = path_as_ind[i]
        d = dist_mat[c1][c2]
        cum = d + prevT
        time.append(cum)
        prevT = cum

    return time


def get_distance_matrix(dataF):
    nCol = len(dataF.columns)
    dist_mat = np.full((nCol, nCol), math.inf)

    for i in range(0, nCol - 1):
        v1 = dataF.iloc[:, i].to_numpy()

        for j in range(i + 1, nCol):
            v2 = dataF.iloc[:, j].to_numpy()
            dif = v1 - v2
            d = np.dot(dif, dif)
            dist_mat[i][j] = d
            dist_mat[j][i] = d

    return dist_mat


def get_cosine_distance_matrix(dataF):
    nCol = len(dataF.columns)
    dist_mat = np.full((nCol, nCol), math.inf)

    for i in range(0, nCol - 1):
        v1 = dataF.iloc[:, i].to_numpy()

        for j in range(i + 1, nCol):
            v2 = dataF.iloc[:, j].to_numpy()
            d = distance.cosine(v1, v2)

            # d = d * d
            dist_mat[i][j] = d
            dist_mat[j][i] = d

    return dist_mat


def write_distance_matrix(input_file, output_file):
    dataF = pd.read_table(input_file, index_col=0)
    dist = get_distance_matrix(dataF)
    cols = dataF.columns.values.tolist()
    distF = pd.DataFrame(data=dist, index=cols, columns=cols)
    distF.to_csv(output_file, sep="\t")


def write_cosine_distance_matrix(input_file, output_file):
    dataF = pd.read_table(input_file, index_col=0)
    dist = get_cosine_distance_matrix(dataF)
    cols = dataF.columns.values.tolist()
    distF = pd.DataFrame(data=dist, index=cols, columns=cols)
    distF.to_csv(output_file, sep="\t")


def get_traj_and_time(dist_file: str, source_cell: str, target_cell: str):
    distF = pd.read_table(dist_file, index_col=0)
    return get_traj_and_time_from_dataframe(distF)


def get_traj_and_time_from_dataframe(distF, source_cell: str, target_cell: str):
    cols = distF.columns.values.tolist()
    dist_mat = distF.to_numpy()
    src_ind = cols.index(source_cell)
    trg_ind = cols.index(target_cell)

    prev = get_shortest_path(dist_mat, src_ind, trg_ind)
    path_as_ind = get_path_as_indices(prev, src_ind, trg_ind)
    time = get_pseudotime(path_as_ind, dist_mat)
    path_as_cells = get_path_as_cell_names(path_as_ind, cols)
    return path_as_cells, time


def get_consecutive_distances(distF, traj):
    indices = get_indices(distF.columns.values.tolist(), traj)
    cons = []
    mat = distF.to_numpy()
    for i in range(len(traj) - 1):
        cons.append(mat[indices[i], indices[i+1]])
    return cons


def get_indices(cells, traj):
    indices = []
    for i in range(0, len(traj)):
        indices.append(cells.index(traj[i]))
    return indices


def get_traj_sorted_values(all_vals, traj_ind):
    vals = []

    for i in traj_ind:
        vals.append(all_vals[i])

    return vals


def write_changes(input_path, traj, time, output_path):
    dataF = pd.read_table(input_path, index_col=0)
    cols = dataF.columns.values.tolist()
    genes = dataF.index.values.tolist()
    indices = get_indices(cols, traj)
    mat = dataF.to_numpy()
    pvals = {}

    for i in range(0, len(genes)):
        gene = genes[i]
        vals = get_traj_sorted_values(mat[i], indices)
        tau, p_value = stats.kendalltau(vals, time)

        if not math.isnan(p_value):
            if tau < 0:
                p_value = -p_value

            pvals[gene] = p_value

    with open(output_path, 'w') as f:
        f.write("Gene\tSignedP")
        for key, value in pvals.items():
            f.write('\n%s\t%e' % (key, value))


def k_rank_svd(input_path: str, k: int, output_path: str):
    dataF = pd.read_table(input_path, index_col=0)
    cells = dataF.columns.values.tolist()
    genes = dataF.index.values.tolist()
    mat = dataF.to_numpy()

    u, s, vt = randomized_svd(mat, n_components=k, random_state=None)
    mat_reduced = np.matmul(np.matmul(u, np.diag(s)), vt)

    # Set negative values to 0
    for i in range(len(genes)):
        for j in range(len(cells)):
            if mat_reduced[i][j] < 0:
                mat_reduced[i][j] = 0

    dataF_red = pd.DataFrame(mat_reduced, columns=cells, index=genes)
    dataF_red.to_csv(output_path, sep="\t")


def transpose_table(input_path: str, output_path: str):
    dataF = pd.read_table(input_path, index_col=0)
    dataF.transpose().to_csv(output_path, sep="\t")


def get_cell_closest_to_the_average_loc(input_file, subset):
    dataF = pd.read_table(input_file, index_col=0)
    mid = np.zeros(len(dataF))
    for cell in subset:
        mid += dataF[cell].to_numpy()
    mid /= len(subset)

    closest = None
    min_dist = math.inf

    for cell in subset:
        vec = dataF[cell].to_numpy()
        dist = distance.cosine(mid, vec)
        if dist < min_dist:
            min_dist = dist
            closest = cell

    return closest