
from sklearn.decomposition import NMF
import scipy.stats
import pandas as pd
import numpy as np
import networkx as nx
import os

from ppi_tools import build_PPI

def main(output_folder):
    os.makedirs(output_folder,exist_ok=True)
    
    ppi = build_PPI()

    nodelist = list(ppi.nodes())
    adjacency = nx.adjacency_matrix(ppi, nodelist=nodelist)

    degree = np.array([ppi.degree(i) for i in nodelist]).reshape((-1,1))

    evals,evecs = spectrum(adjacency, 50)
    sle_factors = pd.DataFrame(evecs/evecs.std(axis=0),index=nodelist,columns=[f"SN_Factor_{i}" for i in range(evecs.shape[1])])
    # random walk (column-normalized adjacency) eigenvectors are simply a degree-scaling of the symmetric Laplacian eigenvectors
    rwd_factors = pd.DataFrame((degree*evecs)/(degree*evecs).std(axis=0),index=nodelist,columns=[f"RW_Factor_{i}" for i in range(evecs.shape[1])])

    nmf_factors = make_nmf_factors(adjacency,nodelist)
    
    rwd_factors.to_csv(f"{output_folder}/rwd_factors.csv")
    uniformize(rwd_factors).to_csv(f"{output_folder}/uniform_rwd_factors.csv")
    normalize(rwd_factors).to_csv(f"{output_folder}/normal_rwd_factors.csv")

    sle_factors.to_csv(f"{output_folder}/sle_factors.csv")
    uniformize(sle_factors).to_csv(f"{output_folder}/uniform_sle_factors.csv")
    normalize(sle_factors).to_csv(f"{output_folder}/normal_sle_factors.csv")
    
    nmf_factors.to_csv(f"{output_folder}/nmf_factors.csv")


def make_nmf_factors(adjacency,nodelist):
    model = NMF(n_components=50, init='random', random_state=0,max_iter=10000, verbose=True)
    W = model.fit_transform(adjacency)
    nmf_factors = pd.DataFrame(W,
                            index=nodelist,
                            columns=[f"NMF_Factor_{i}" for i in range(W.shape[1])]
                            )
    
    return nmf_factors

def spectrum(adjacency, n_evectors):
    """
    Compute the top eigenvalues and eigenvectors of the normalized adjacency matrix
    derived from a given graph's adjacency matrix, using spectral graph theory.

    This function performs the following steps:
    1. Constructs the normalized adjacency matrix of the graph.
    2. Computes the eigenvalues and eigenvectors of the normalized adjacency matrix
       using `scipy.sparse.linalg.eigsh`.
    3. Returns the eigenvalues and eigenvectors sorted by descending eigenvalue order.

    Parameters
    ----------
    adjacency : np.ndarray or scipy.sparse matrix
        The adjacency matrix of the graph, where `adjacency[i, j]` represents
        the edge weight between node `i` and node `j`. It can be a dense NumPy array
        or a sparse matrix from `scipy.sparse`.
    n_evectors : int
        The number of eigenvectors and eigenvalues to compute (i.e., the top `n_evectors`).

    Returns
    -------
    e : np.ndarray
        The top `n_evectors` eigenvalues of the normalized adjacency matrix, sorted in descending order.
    evecs : np.ndarray
        The corresponding eigenvectors, where each column `evecs[:, i]` is the eigenvector
        associated with eigenvalue `e[i]`, sorted in descending order of eigenvalue.

    Notes
    -----
    - The normalized adjacency matrix is defined as `D^(-1/2) * A * D^(-1/2)`, where `A`
      is the original adjacency matrix and `D` is the diagonal degree matrix with
      `D[i, i] = sum(A[i, :])`. This is derived from the graph's Laplacian matrix.
    - The normalized adjacency is related to the Laplacian matrix:
      `I - L_normalized`, where `L_normalized` is the normalized Laplacian.
    - The function computes the eigenvectors of the normalized adjacency matrix,
      which are often used in spectral clustering or other graph-based algorithms.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csr_matrix
    >>> adjacency = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    >>> spectrum(adjacency, 2)
    (array([1.618, 0.618]), array([[0.5, -0.5], [0.707, 0.707], [-0.5, -0.5]]))
    """

    # Number of nodes in the graph (size of adjacency matrix)
    n_nodes = adjacency.shape[0]

    # Compute the degree of each node (sum of each row in the adjacency matrix)
    d = adjacency @ np.ones(n_nodes)

    # Construct the diagonal degree matrix D
    D = scipy.sparse.diags(d)

    # Compute D^(-1/2), the inverse square root of the degree matrix
    D_neg_half = scipy.sparse.diags(d**(-1/2))

    # Construct the normalized adjacency matrix: D^(-1/2) * A * D^(-1/2)
    normalized_adjacency = D_neg_half @ adjacency @ D_neg_half

    # Calculate the eigenvectors and eigenvalues of the normalized adjacency matrix
    e, evecs = scipy.sparse.linalg.eigsh(normalized_adjacency, k=n_evectors)

    # The eigenvalues from eigsh are returned in ascending order,
    # so we reverse them to get descending order (high eigenvalue = low Laplacian eigenvalue)
    e = e[::-1]
    evecs = evecs[:, ::-1]

    return e, evecs


def uniformize(factors):
    factors_uniform = scipy.stats.rankdata(factors,axis=0)/factors.shape[0]
    factors_uniform = pd.DataFrame(factors_uniform,index=factors.index,columns=[f"Uniform_{factor_name}" for factor_name in factors.columns])
    return factors_uniform


def normalize(factors):
    factors_uniform = uniformize(factors)
    
    eps = 0.001
    
    # we move the uniform ever so slightly 
    # away from 0 to 1 as this avoids creating
    # infinite values at the extreme points
    factors_normal = scipy.stats.norm.isf(1-((1-eps)*factors_uniform+eps/2))
    factors_normal = pd.DataFrame(factors_normal,index=factors.index,columns=[f"Normal_{factor_name}" for factor_name in factors.columns])
    
    return factors_normal

if __name__ == "__main__":
    
    main("data/intermediate/factor_dataframes")