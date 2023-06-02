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

# Generates transition probabilities from the similarities
def get_probability_matrix(sim_mat):
    sums = np.sum(sim_mat, axis=1)
    p_mat = (sim_mat.T / sums).T
    return p_mat

# Generates the flows between pairs of cells using transition probabilities
# and the current heat distribution
def get_flow_matrix(p_matrix, heat):
    return (heat * p_matrix.T).T

# Next round of heat distribution is determined by the incoming fow to each cell
def get_next_heat(flow):
    return np.sum(flow, axis=0)

def generate_flow_matrices(p_matrix, source_indicator, target_indicator, k:int):
    # Heat flow at each step
    flow_list = []

    # Heat that arrives to the target(s) at each step. That heat is absorbed by target(s),
    # hence we say the heat sinks.
    sink_list = []

    # Heat of cells
    heat = np.zeros(len(source_indicator))

    # Source cells are hot (1), all other cells are cold (0).
    heat[source_indicator] = 1

    # The differential sinking heat in last step. If the parameter k (total number of steps) is given as -1, then the simulation goes on until this differential becomes negative.
    diff = 1

    while (k < 0 < diff) or (k > len(flow_list)):
        # Get the next flows using the probability matrix and current heat distribution
        flow = get_flow_matrix(p_matrix, heat)
        flow_list.append(flow)

        # Calculate the new heat distribution based on the new flows
        heat = get_next_heat(flow)

        # Record how much of this flow reached target(s), and sank
        sink_list.append(np.sum(heat[target_indicator]))

        # Check the trend in sinking heat
        if len(sink_list) > 1:
            diff = sink_list[-1] - sink_list[-2]

    return flow_list, sink_list

# For the heat that has arrived at iteration k, this function generates the list of heat distribution at each step. It also returns a cumulative flow matrix that focuses on the heat that arrived to target.
def generate_weights(flow_matrices, target_indicator, arrive_index:int):
    if arrive_index == -1:
        arrive_index = len(flow_matrices) - 1

    # Part of the flow that arrives to target
    sub_flow = np.zeros(flow_matrices[0].shape)

    # Number of cells
    cell_size = len(target_indicator)

    # The weights matrix where each row is for one turn of heat distribution
    # Weight: The heat that will arrive to the target(s)
    weights_mat = np.zeros((arrive_index + 2, cell_size))

    # The ratios array indicates how much of the in-flow will arrive to the target(s)
    ratios = np.zeros(cell_size)

    # The last step in-flow to the targets will arrive to targets in full, so their ratios are 1
    ratios[target_indicator] = 1

    # The array of in-flows to each cell
    inflow_sums = np.sum(flow_matrices[arrive_index], axis = 0)

    # Weights at the last step is in-flows to the targets and zero for everything else
    weights = inflow_sums * ratios

    # Put this weights array to the last row of weights matrix
    weights_mat[arrive_index + 1, :] = weights

    # Now we will calculate the rest of the weight arrays by going backwards in time
    index = arrive_index

    while index >= 0:
        # Calculate prior weights
        weights = np.matmul(flow_matrices[index], ratios.T)

        # Identify the components of flows at this step that arrived to the target(s) at arrive_index
        # Add those flows to the cumulative matrix
        sub_flow += flow_matrices[index] * ratios

        # Find out-flows for each cell
        outflow_sums = np.sum(flow_matrices[index], axis = 1)

        # Calculate the new ratios of emitted heat that arrive at the target(s)
        ratios = weights / outflow_sums

        # Targets have zero outflows. That generates nans in the previous division. We replace them with zeros here.
        ratios[np.isnan(ratios)] = 0

        # Save new weights
        weights_mat[index, :] = weights

        # Go one back in time
        index = index - 1

    return weights_mat, sub_flow

# This is the cumulative version of the above function. It is for the heat that arrived at iteration k or earlier.
def generate_cumulative_weights(flow_matrices, target_indicator, max_arrive_index):

    # If arrive_index is negative, use the last flow matrix
    if max_arrive_index < 0:
        max_arrive_index = len(flow_matrices) - 1

    # This is the double cumulative sub-flow matrix
    sub_flow_cum = np.zeros(flow_matrices[0].shape)

    # We will generate weight arrays for the heat that arrives at arrive index, then for the heat that arrives at arrive_index-1 and so on.
    index = max_arrive_index
    m_list = []
    while index >= 0:
        # Get the weight arrays and sub-flow matrix for the heat that arrives to target(s) at step "index"
        w, f = generate_weights(flow_matrices, target_indicator, index)

        # Save the weight matrix
        m_list.append(w)

        # Add the cumulative sub-flows to the double-cumulative sub-flow matrix
        sub_flow_cum += f

        # Go back in arrival time
        index -= 1

    # Initialize cumulative weights
    weights = np.zeros(m_list[0].shape)

    # Integrate the weight matrices. After integration, weights are the heat on these cells that will arrive at the target(s) at step arrive_index or earlier
    for i in range(len(weights)):
        for j in range(min(max_arrive_index - i + 2, max_arrive_index + 1)):
            weights[i] += m_list[j][i]

    return sub_flow_cum
