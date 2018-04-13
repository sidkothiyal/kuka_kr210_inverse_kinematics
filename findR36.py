from sympy import Matrix, symbols, sin, cos, atan2, pi, simplify


def rot_mat(a, q):
	T = Matrix([[	cos(q),       -sin(q),	     0],
		[sin(q)*cos(a), cos(q)*cos(a), -sin(a)],
		[sin(q)*sin(a), cos(q)*sin(a),  cos(a)]])
	return T

if __name__ == "__main__":
	q4, q5, q6 = symbols("q4:7")
	a3, a4, a5 = symbols("a3:6")

	R3_6 = simplify(rot_mat(a3, q4) * rot_mat(a4, q5) * rot_mat(a5, q6)).subs({a3: -pi/2, a4: pi/2, a5: -pi/2})
	print R3_6
