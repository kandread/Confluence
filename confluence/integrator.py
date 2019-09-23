"""

Algorithms for integrating SWOT-estimated discharge.

"""

from functools import partial
from scipy.optimize import minimize
import numpy as np
import pandas as pd


class MeanFlow:
    """Mean flow discharge integrator."""

    def __init__(self, n, A0, H, W, S, A, dA, routing_table):
        self.nreaches = len(A0)
        self.A0 = np.array(A0)
        self.n = np.array(n)
        self.data = (H, W, S, A, dA)
        self.rivs = self._riverTopology(routing_table)

    def _riverTopology(self, routing_table):
        """Create river topology from routing table. River topology is
        a dictionary with each entry being a list of upstream reaches for
        the specific reach."""
        rout = pd.read_csv(routing_table)
        rout = rout.dropna(axis=0)[['GridID', 'HydroID', 'NextDownID']]
        rivs = {k: [] for k in rout['GridID']}
        hyd = {row[1]['HydroID']: row[1]['GridID'] for row in rout.iterrows()}
        for row in rout.iterrows():
            if row[1]['NextDownID'] > 0:
                rivs[hyd[row[1]['NextDownID']]].append(hyd[row[1]['HydroID']])
        idx = {k: i for i, k in enumerate(rivs)}
        rivs = {i: [idx[l] for l in rivs[k]] for i, k in enumerate(rivs)}
        return rivs

    def objective(self, x):
        """Objective function for mean-flow integration."""
        n = np.array(x[1::2])
        A0 = np.array(x[0::2])
        H, W, S, A, dA = self.data
        w = np.mean(W, axis=0)
        h = np.mean(H, axis=0)
        if dA is None:
            dA = np.array([(w[r] + W[np.argmin(A[:, r]), r]) / 2 * (h[r] - H[np.argmin(A[:, r]), r])
                           for r in range(self.nreaches)]).T
        else:
            dA = np.mean(dA, axis=0)
        Q0 = 1 / self.n * (self.A0 + dA)**(5 / 3) * w**(-2 / 3) * np.mean(
            S, axis=0)**(1 / 2)
        Q = 1 / n * (A0 + dA)**(5 / 3) * w**(-2 / 3) * np.mean(S,
                                                               axis=0)**(1 / 2)
        return np.sum((Q - np.mean(Q0))**2)

    def constraint(self, i, x):
        """Constrain discharge to increase downstream."""
        n = np.array(x[1::2])
        A0 = np.array(x[0::2])
        H, W, S, A, dA = self.data
        w = np.mean(W, axis=0)
        h = np.mean(H, axis=0)
        if dA is None:
            dA = np.array([(w[r] + W[np.argmin(A[:, r]), r]) / 2 * (h[r] - H[np.argmin(A[:, r]), r])
                           for r in range(self.nreaches)]).T
        else:
            dA = np.mean(dA, axis=0)
        Q = 1 / n * (A0 + dA)**(5 / 3) * w**(-2 / 3) * np.mean(S,
                                                               axis=0)**(1 / 2)
        return Q[i] - Q[i + 1]

    def integrate(self):
        """Integrate discharge by forcing mean annual
        flow to increase downstream."""
        lbnds = [0.0] * (self.nreaches * 2)
        ubnds = [0.0] * (self.nreaches * 2)
        lbnds[1::2] = [0.01] * self.nreaches
        ubnds[1::2] = [0.09] * self.nreaches
        ubnds[::2] = [1e4] * self.nreaches
        bnds = list(zip(lbnds, ubnds))
        x0 = [0.0] * (self.nreaches * 2)
        x0[1::2] = self.n
        x0[::2] = self.A0
        cons = [{
            'type': 'ineq',
            'fun': partial(self.constraint, i)
        } for i in range(self.nreaches - 1)]
        solution = minimize(self.objective,
                            x0,
                            method='trust-constr',
                            bounds=bnds,
                            tol=1e-3,
                            constraints=cons)
        return solution.x[::2], solution.x[1::2]
