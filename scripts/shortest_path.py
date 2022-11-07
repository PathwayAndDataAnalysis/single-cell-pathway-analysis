import pandas as pd
import numpy as np
import math


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


def write_distance_matrix(input_file, output_file):
    dataF = pd.read_table(input_file, index_col=0)
    dist = get_distance_matrix(dataF)
    cols = dataF.columns.values.tolist()
    distF = pd.DataFrame(data=dist, index = cols, columns = cols)
    distF.to_csv(output_file, sep="\t")


def get_traj_and_time(dist_file:str, source_cell:str, target_cell:str):
    distF = pd.read_table(dist_file, index_col=0)
    cols = distF.columns.values.tolist()
    dist_mat = distF.to_numpy()
    src_ind = cols.index(source_cell)
    trg_ind = cols.index(target_cell)

    prev = get_shortest_path(dist_mat, src_ind, trg_ind)
    path_as_ind = get_path_as_indices(prev, src_ind, trg_ind)
    time = get_pseudotime(path_as_ind, dist_mat)
    path_as_cells = get_path_as_cell_names(path_as_ind, cols)
    return path_as_cells, time


