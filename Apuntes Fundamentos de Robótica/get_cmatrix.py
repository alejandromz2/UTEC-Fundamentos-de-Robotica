def get_cmatrix(M,q_v,dq_v):

	n = len(q_v)

	if n == 2:

		q1 = q_v[0]; q2 = q_v[1];
		dq1 = dq_v[0]; dq2 = dq_v[1];

		c111 = (sp.diff(M[0,0],q1)+sp.diff(M[0,0],q1)-sp.diff(M[0,0],q1))/2
		c112 = (sp.diff(M[0,0],q2)+sp.diff(M[0,1],q1)-sp.diff(M[0,0],q1))/2
		c121 = c112
		c122 = (sp.diff(M[0,1],q2)+sp.diff(M[0,1],q2)-sp.diff(M[1,1],q1))/2
		c211 = (sp.diff(M[1,0],q1)+sp.diff(M[1,0],q1)-sp.diff(M[0,0],q2))/2
		c212 = (sp.diff(M[1,0],q2)+sp.diff(M[1,1],q1)-sp.diff(M[0,1],q2))/2
		c221 = c212
		c222 = (sp.diff(M[1,1],q1)+sp.diff(M[1,1],q2)-sp.diff(M[1,1],q2))/2

		c11 = c111*dq1 + c112*dq2
		c12 = c121*dq1 + c122*dq2
		c21 = c211*dq1 + c212*dq2
		c22 = c221*dq1 + c222*dq2

		C = sp.Matrix([[c11,c12],[c21,c22]])

		return C

	else:

		return None