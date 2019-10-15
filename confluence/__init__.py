"""

Implementation of Confluence framework for integrating and diagnosing
SWOT-estimated discharge.

"""

import numpy as np


def removeFlagged(H, W, S, dA):
    """Clean up data."""
    i = np.where(np.all(H > -1000, axis=1))[0]
    H = H[i, :]
    i = np.where(np.all(W > -1000, axis=1))[0]
    W = W[i, :]
    i = np.where(np.all(S > -1000, axis=1))[0]
    S = S[i, :]
    if dA is not None:
        dA = dA[i, :]
    return H, W, S, dA
