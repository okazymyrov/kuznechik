load ./Kuznechik_2_8.sage
load ./Kuznechik.sage
load ./Kuznechik_table.sage

def generate_mds_matrix(**kwargs):
	r'''
		Generate MDS matrix for the L transformation.

		INPUT::

			- ``method`` -- 'basic', 'mathematical' and 'alternative' approaches to generate the matrix.
							While 'basic' and 'mathematical' are equivalent approaches, 'alternative' is
							equivalent to transposed of 'basic'

			- ``output``  -- boolean (default: False), output the matrix ot stdout
	'''
	method = kwargs.get('method','basic')

	K = GF(2^8,'a',modulus=ZZ['x']("x^8+x^7+x^6+x+1"))

	if method == 'basic':
		P = PolynomialRing(K,16,'x')

		# l routine over P
		def l(state):
			f = P("0")
			for i,x in enumerate([148, 32, 133, 16, 194, 192, 1, 251, 1, 192, 194, 16, 133, 32, 148, 1]):
				f += P("({})*({})".format(K.fetch_int(x),state[i]))

			return f

		# R routine over P
		def R(state):
			return [l(state)] + state[:15]

		# L routine over P
		def L(state):
			for i in xrange(16):
				state = R(state)

			return state

		state = list(P.gens())[::-1]
		f = L(state)
		m,_ = Sequence(f).coefficient_matrix()

		m = m*Permutation([g for g in xrange(1,17)][::-1]).to_matrix()

		MDS = matrix(ZZ,16)

		for i in xrange(16):
			for j in xrange(16):
				MDS[i,j] = m[i,j].integer_representation()

	elif method == 'mathematical':
		t = companion_matrix([K.fetch_int(g) for g in [1,148, 32, 133, 16, 194, 192, 1, 251, 1, 192, 194, 16, 133, 32, 148, 1][::-1]], format='top')^16

		MDS = matrix(ZZ,16)

		for i in xrange(16):
			for j in xrange(16):
				MDS[i,j] = t[i,j].integer_representation()

	elif method == 'alternative':
		t = companion_matrix([K.fetch_int(g) for g in [1,148, 32, 133, 16, 194, 192, 1, 251, 1, 192, 194, 16, 133, 32, 148, 1][::-1]], format='left')^16

		MDS = matrix(ZZ,16)

		for i in xrange(16):
			for j in xrange(16):
				MDS[i,j] = t[i,j].integer_representation()

	else:
		raise ValueError("Unknown 'method'")

	if kwargs.get('output',False):
		sys.stdout.write("MDS:\n[\n")
		for i in xrange(16):
			sys.stdout.write("\t[")
			for j in xrange(16):
				if j == 15:
					sys.stdout.write("0x{0:02X}],\n".format(MDS[i,j]))
				else:
					sys.stdout.write("0x{0:02X},".format(MDS[i,j]))
		sys.stdout.write("]\n")

	return MDS

def generate_inverse_mds_matrix(**kwargs):
	r'''
		Generate MDS matrix for the L transformation.

		INPUT::

			- ``method`` -- 'basic', 'mathematical' and 'alternative' approaches to generate the matrix.
							While 'basic' and 'mathematical' are equivalent approaches, 'alternative' is
							equivalent to transposed of 'basic'

			- ``output``  -- boolean (default: False), output the matrix ot stdout
	'''

	K = GF(2^8,'a',modulus=ZZ['x']("x^8+x^7+x^6+x+1"))

	MDS = generate_mds_matrix(method=kwargs.get('method','basic'))

	T = matrix(K,16)

	for i in xrange(16):
		for j in xrange(16):
			T[i,j] = K.fetch_int(MDS[i,j])

	T = T.inverse()

	MDS_inv = matrix(ZZ,16)

	for i in xrange(16):
		for j in xrange(16):
			MDS_inv[i,j] = T[i,j].integer_representation()

	if kwargs.get('output',False):
		sys.stdout.write("MDS_inv:\n[\n")
		for i in xrange(16):
			sys.stdout.write("\t[")
			for j in xrange(16):
				if j == 15:
					sys.stdout.write("0x{0:02X}],\n".format(MDS_inv[i,j]))
				else:
					sys.stdout.write("0x{0:02X},".format(MDS_inv[i,j]))
		sys.stdout.write("]\n")

	return MDS_inv

def generate_tables(**kwargs):
	r'''
		Generate MDS matrix for the L transformation.

		INPUT::
			- ``filename``  -- filename for the output file (default: stdout)
	'''

	K = GF(2^8,'a',modulus=ZZ['x']('x^8+x^7+x^6+x+1'))
	G = matrix(K,M.nrows())

	for i in xrange(M.nrows()):
		for j in xrange(M.ncols()):
			G[i,j] = K.fetch_int(M[i,j])

	G_inv = matrix(K,M_inv.nrows())

	for i in xrange(M_inv.nrows()):
		for j in xrange(M_inv.ncols()):
			G_inv[i,j] = K.fetch_int(M_inv[i,j])

	def from_state(state):
		return ZZ([g for g in state[::-1]],2^8)

	def L(state):
		return [g.integer_representation() for g in (G*matrix(K,16,1,[K.fetch_int(g) for g in state])).list()]

	def L_inv(state):
		return [g.integer_representation() for g in (G_inv*matrix(K,16,1,[K.fetch_int(g) for g in state])).list()]

	def to_state(state):
		return [g for g in ZZ(state).digits(base=2^8,padto=16)[::-1]]

	filename = kwargs.get('filename',None)

	if filename is None:
		output = sys.stdout
	else:
		output = open(filename, "w")

	LS = [[] for i in xrange(16)]
	for i in xrange(16):
		for j in xrange(256):
			state = to_state(Sbox[j]<<(8*(15-i)))
			LS[i].append(from_state(L(state)))

	LS_inv = [[] for i in xrange(16)]
	for i in xrange(16):
		for j in xrange(256):
			state = to_state(Sbox_inv[j]<<(8*(15-i)))
			LS_inv[i].append(from_state(L_inv(state)))

	output.write("# Tables for the LS transformation\nLS = [\n")
	for i in xrange(16):
		output.write("[\n\t")
		for j,s in enumerate(LS[i]):
			if (j+1) % 4 == 0:
				if j == 255:
					output.write("0x{0:032X}\n".format(s))
				else:
					output.write("0x{0:032X},\n\t".format(s))
			else:
				output.write("0x{0:032X},".format(s))
		if i == 15:
			output.write("]\n")
		else:
			output.write("],\n")
	output.write("]\n")

	output.write("\n")

	output.write("# Tables for the inverse LS transformation\nLS_inv = [\n")
	for i in xrange(16):
		output.write("[\n\t")
		for j,s in enumerate(LS_inv[i]):
			if (j+1) % 4 == 0:
				if j == 255:
					output.write("0x{0:032X}\n".format(s))
				else:
					output.write("0x{0:032X},\n\t".format(s))
			else:
				output.write("0x{0:032X},".format(s))
		if i == 15:
			output.write("]\n")
		else:
			output.write("],\n")
	output.write("]\n")

	if filename is not None:
		output.close()

def test_matrix_generation():
	r'''
		Test method for generation MDS matricies
	'''
	MDS1 = generate_mds_matrix(method='basic')
	MDS2 = generate_mds_matrix(method='mathematical')
	MDS3 = generate_mds_matrix(method='alternative')

	print "{0} {1}".format("".center(14,'-'),"".center(6,'-'))
	print "{0} {1}".format("method".center(14),"matrix".center(6))
	print "{0} {1}".format("".center(14,'-'),"".center(6,'-'))
	print "{0} {1}".format("basic".center(14),"MDS1".center(6))
	print "{0} {1}".format("mathematical".center(14),"MDS2".center(6))
	print "{0} {1}".format("alternative".center(14),"MDS3".center(6))
	print "{0} {1}".format("".center(14,'-'),"".center(6,'-'))

	print "MDS1 == MDS2\t: {0}".format(MDS1==MDS2)
	print "MDS1 == MDS3^t\t: {0}".format(MDS1==MDS3.transpose())
	print "MDS2 == MDS3^t\t: {0}".format(MDS2==MDS3.transpose())

def test_performance():
	r'''
		Test performance of the proposed implementations
	'''
	K = ZZ.random_element(2^256)
	M = ZZ.random_element(2^128)

	space = 1<<5

	cipher = Kuznechik_2_8(key=K)

	t1 = cputime()

	for _ in xrange(space):
		M = cipher.encrypt(M)

	t2 = cputime()

	print "Kuznechik_2_8\t: {0}".format(t2-t1)

	cipher = Kuznechik(key=K)

	t1 = cputime()

	for _ in xrange(space):
		M = cipher.encrypt(M)

	t2 = cputime()

	print "Kuznechik\t: {0}".format(t2-t1)

	cipher = Kuznechik_table(key=K)

	t1 = cputime()

	for _ in xrange(space):
		M = cipher.encrypt(M)

	t2 = cputime()

	print "Kuznechik_table\t: {0}".format(t2-t1)
