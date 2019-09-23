"""

Implementation of Confluence framework for integrating and diagnosing
SWOT-estimated discharge.

"""

import numpy as np


def removeFlagged(H, W, S):
    """Clean up data."""
    i = np.where(np.all(H > -1000, axis=1))
    H = H[i, :]
    i = np.where(np.all(W > -1000, axis=1))
    W = W[i, :]
    i = np.where(np.all(S > -1000, axis=1))
    S = S[i, :]
    return H, W, S


