r"""
    Implementaion of Kuznechik over GF(2^8)

AUTHORS:

- Oleksandr Kazymyrov (2014-08-08): initial version

"""

#*****************************************************************************
#       Copyright (C) 2014 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***************************************************************************** 

load ./Data.sage
load ./Tables.sage

class Kuznechik_table(SageObject):
	def __init__(self, **kwargs):
		r'''
			Implementation of the block cipher Kuznechik over GF(2^8)

			INPUT::

				- ``key`` 		-- the ecnryption key
		'''
		
		self._kl 	  = 256 # bits
		self._bl 	  = 128  # bits
		self._rks 	  = [0 for _ in xrange(10)]
		self._rks_inv = [0 for _ in xrange(10)]
	
		self._key 	= kwargs.get('key',None)

		self._Sbox = Sbox
		self._Sbox_inv = Sbox_inv

		if self._key is not None:
			self.__KeyExpansion(self._key)
		else:
			raise ValueError("key must be presented")

	def decrypt(self,ciphertext):
		r'''
			Decryption routine
		'''
		state = self.__to_state(ciphertext)

		state = self.__S(state)

		for i in xrange(9,0,-1):
			state = self.__LS_inv(state)
			state = self.__X(self._rks_inv[i],state)

		state = self.__S_inv(state)
		state = self.__X(self._rks_inv[0],state)

		return self.__from_state(state)

	def encrypt(self,plaintext):
		r'''
			Encryption routine
		'''

		state = self.__to_state(plaintext)

		state = self.__X(self._rks[0],state)

		for i in xrange(1,10):
			state = self.__LS(state)
			state = self.__X(self._rks[i],state)

		return self.__from_state(state)

	def __from_state(self,state):
		r'''
			Convert the iternal representation to the corresponding integer 
		'''	
		return ZZ([g for g in state[::-1]],2^8)

	def __KeyExpansion(self,K):
		r'''
			Generate subkeys for all round of encryption based on the given master key
		'''
		C = [self.__LS(self.__S_inv(self.__to_state(i))) for i in xrange(1,33)]

		self._rks[0] = self.__to_state(ZZ(K).digits(2^128,padto=2)[1])
		self._rks[1] = self.__to_state(ZZ(K).digits(2^128,padto=2)[0])

		l = self._rks[0]
		r = self._rks[1]

		for i in xrange(4):
			for j in xrange(0,8,2):
				r = self.__X(r,self.__LS(self.__X(C[i*8+j],l)))
				l = self.__X(l,self.__LS(self.__X(C[i*8+j+1],r)))
			self._rks[(i+1)*2]   = l
			self._rks[(i+1)*2+1] = r

		self._rks_inv = [self._rks[0]] + [self.__LS_inv(self.__S(k)) for k in self._rks[1:]]

	def __LS(self,state):
		r'''
			The LS routine
		'''
		t = 0

		for i in xrange(len(state)):
			t = t ^^ LS[i][state[i]]

		return self.__to_state(t)

	def __LS_inv(self,state):
		r'''
			The inverse LS routine
		'''
		t = 0

		for i in xrange(len(state)):
			t = t ^^ LS_inv[i][state[i]]

		return self.__to_state(t)

	def __S(self,state):
		r'''
			The S routine
		'''	
		return [self._Sbox[g] for g in state]

	def __S_inv(self,state):
		r'''
			The S^{-1} routine
		'''	
		return [self._Sbox_inv[g] for g in state]

	def selfTesting(self):
		r'''
			Testing the class using data from the specification
		'''

		test_K = 0x8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef

		test_rks = [0x8899aabbccddeeff0011223344556677, 0xfedcba98765432100123456789abcdef, 0xdb31485315694343228d6aef8cc78c44, 0x3d4553d8e9cfec6815ebadc40a9ffd04, 
			 		0x57646468c44a5e28d3e59246f429f1ac, 0xbd079435165c6432b532e82834da581b, 0x51e640757e8745de705727265a0098b1, 0x5a7925017b9fdd3ed72a91a22286f984, 
			 		0xbb44e25378c73123a5f32f73cdb6e517, 0x72e9dd7416bcf45b755dbaa88e4a4043]

		rks     = self._rks[:]
		rks_inv = self._rks_inv[:]
		self.__KeyExpansion(test_K)

		if ([self.__from_state(k) for k in self._rks] != test_rks):
			raise TypeError("Self testing of KeyExpansion fail")

		M = 0x1122334455667700ffeeddccbbaa9988
		C = 0x7f679d90bebc24305a468d42b9d4edcd

		if (self.encrypt(M) != C):
			raise TypeError("Self testing of encryption fail")

		if (self.decrypt(C) != M):
			raise TypeError("Self testing of decryption fail")		

		self._rks 	  = rks[:]
		self._rks_inv = rks_inv[:]

	def __to_state(self,state):
		r'''
			Convert an integer to the iternal representation
		'''	
		return [g for g in ZZ(state).digits(base=2^8,padto=16)[::-1]]

	def __X(self,key,state):
		r'''
			The X routine
		'''
		return [g[0]^^g[1] for g in zip(key,state)]
