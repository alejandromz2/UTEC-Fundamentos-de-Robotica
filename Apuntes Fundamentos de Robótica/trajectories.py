import numpy as np
import sympy as sp

def get_coeffs(ec_m, v_m, type = 's'):

	nv = v_m.shape[0]
	n_ec = ec_m.shape[0]

	coeffs = sp.zeros(n_ec,nv)

	for i in range(n_ec):
		for j in range(nv):
			coeffs[i,j] = (ec_m[i].evalf()).coeff(v_m[j,0])

	if type == 'n':
		n_coeffs = np.array(coeffs.tolist()).astype(np.float64)
		return n_coeffs
	elif type == 's':
		return coeffs