import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.spatial import distance


# Calculates the similarity matrix based on distances
def get_similarity_matrix(m:np.array, sigmasq:float):
    nRow, nCol = m.shape
    sim_mat = np.full((nCol, nCol), 0.0)

    for i in range(0, nCol):
        for j in range(i, nCol):

            # Squared Euclidean distance
            dif = m[:, i] - m[:, j]
            d = np.dot(dif, dif)

            # Cosine distance
            # d = distance.cosine(m[:, i], m[:, j])

            # Use Gaussian similarity kernel to convert distance into similarity
            s = np.exp(-d / sigmasq)

            # Save the similarity on the result matrix, which is symmetric
            sim_mat[i][j] = s
            sim_mat[j][i] = s

    return sim_mat
