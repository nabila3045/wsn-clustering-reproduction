"""Regression test: hand-checkable multi-hop example from MODEL_SPEC.md (no setup, 4000-bit payload).
M --20 m--> CH A --40 m--> CH B --60 m--> sink. Expected debits (mJ): M 0.216, A 0.504, B 0.584, total 1.304."""
import sys, numpy as np
sys.path.insert(0, '.')
import simulator as s
c = {**s.CFG, 'initial_energy_j': 0.5}
xy = np.array([[0., 0.], [20., 0.], [60., 0.]]); sink = np.array([120., 0.])      # M, A, B
dist = np.linalg.norm(xy[:, None]-xy[None], axis=2); ds = np.linalg.norm(xy-sink, axis=1); e = np.full(3, 0.5)
debit, parent, role, txb, rxb, dab, ce = s.costs(np.array([1, 2]), np.ones(3, bool), dist, ds, e, c, True, False)
mj = debit*1e3; print('debits (mJ):', mj.round(4), 'total', mj.sum().round(4))
assert np.allclose(mj, [0.216, 0.504, 0.584], atol=1e-9) and abs(mj.sum()-1.304) < 1e-9
assert list(parent) == [1, 2, -1] and list(role) == [1, 2, 3]       # M->A, A->B, B->BS ; B is the relay CH (role 3)
print('PASS')
