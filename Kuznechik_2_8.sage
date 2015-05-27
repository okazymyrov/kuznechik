r"""
    Implementaion of Kuznechik over GF(2^8)

AUTHORS:

- Oleksandr Kazymyrov (2014-08-06): initial version

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

class Kuznechik_2_8(SageObject):
	def __init__(self, **kwargs):
		r'''
			Implementation of the block cipher Kuznechik over GF(2^8)

			INPUT::

				- ``key`` 		-- the ecnryption key
		'''
		
		self._kl 	= 256 # bits
		self._bl 	= 128  # bits
		self._rks 	= [0 for _ in xrange(10)]
	
		self._key 	= kwargs.get('key',None)

		self._K = GF(2^8,'a',modulus=ZZ['x']('x^8+x^7+x^6+x+1'))
		self._P = PolynomialRing(self._K,'x')

		self._Sbox = self._P("a^30*x^254 + a^132*x^253 + a^61*x^252 + a^104*x^251 + a^32*x^250 + a^198*x^249 + a^38*x^248 + a^132*x^247 + a^16*x^246 + a^88*x^245 \
							 + a^148*x^244 + a^59*x^243 + a^9*x^242 + a^132*x^241 + a^32*x^240 + a^112*x^239 + a^190*x^238 + a^59*x^237 + a^199*x^236 + a^85*x^235 \
							 + a^39*x^234 + a^79*x^233 + a^172*x^232 + a^129*x^231 + a^204*x^230 + a^76*x^229 + a^237*x^228 + a^37*x^227 + a^162*x^226 + a^32*x^225 \
							 + a^100*x^224 + a^162*x^223 + a^100*x^222 + a^43*x^221 + a*x^220 + a^23*x^219 + a^72*x^218 + a^81*x^217 + a^9*x^216 + a^249*x^215 \
							 + a^44*x^214 + a^108*x^213 + a^234*x^212 + a^203*x^211 + a^225*x^210 + a^101*x^209 + a^222*x^208 + a^133*x^207 + a^45*x^206 + a^76*x^205 \
							 + a^186*x^204 + a^32*x^203 + a^169*x^202 + a^83*x^201 + a^163*x^200 + a^172*x^199 + a^153*x^198 + a^107*x^197 + a^49*x^196 + a^66*x^195 \
							 + a^86*x^194 + a^42*x^193 + a^17*x^192 + a^222*x^191 + a^103*x^190 + a^4*x^189 + a^199*x^188 + a^230*x^187 + a^92*x^186 + a^191*x^185 \
							 + a^66*x^184 + a^168*x^183 + a^18*x^182 + a^178*x^181 + a^71*x^180 + a^175*x^179 + a^215*x^178 + a^191*x^177 + a^126*x^176 + a^16*x^175 \
							 + a^120*x^174 + a^144*x^173 + a^93*x^172 + a^116*x^171 + a^158*x^170 + a^9*x^169 + a^180*x^168 + a^228*x^167 + a^173*x^166 + a^64*x^165 \
							 + a^31*x^164 + a^143*x^163 + a^91*x^162 + a^219*x^161 + a^36*x^160 + a^190*x^159 + a^190*x^158 + a^34*x^157 + a^47*x^156 + a^82*x^155 \
							 + a^76*x^154 + a^208*x^153 + a^135*x^152 + a^79*x^151 + a^22*x^150 + a^198*x^149 + a^203*x^148 + a^250*x^147 + a^54*x^146 + a^250*x^145 \
							 + a^217*x^144 + a^33*x^143 + a^75*x^142 + a^225*x^141 + a^174*x^140 + a^52*x^139 + a^155*x^138 + a^38*x^137 + a^42*x^136 + a^201*x^135 \
							 + a^117*x^134 + a^57*x^133 + a^55*x^132 + a^125*x^131 + a^105*x^130 + a^96*x^129 + a^215*x^128 + a^35*x^127 + a^69*x^126 + a^151*x^125 \
							 + a^141*x^124 + a^172*x^123 + a^154*x^122 + a^206*x^121 + a^179*x^120 + a^188*x^119 + a^77*x^118 + a^168*x^117 + a^114*x^116 + a^112*x^115 \
							 + a^38*x^114 + a^198*x^113 + a^211*x^112 + a^210*x^111 + a^231*x^110 + a^34*x^109 + a^229*x^108 + a^113*x^107 + a^39*x^106 + a^203*x^105 \
							 + a^13*x^104 + a^228*x^103 + a^104*x^102 + a^202*x^101 + a^241*x^100 + a^195*x^99 + a^46*x^98 + a^216*x^97 + a^37*x^96 + a^174*x^95 \
							 + a^91*x^94 + a^18*x^93 + a^148*x^92 + a^111*x^91 + a^227*x^90 + a^220*x^89 + a^210*x^88 + a^172*x^87 + a^72*x^86 + a^78*x^85 + a^124*x^84 \
							 + a^92*x^83 + a^237*x^82 + a^171*x^81 + a^223*x^80 + a^113*x^79 + a^124*x^78 + a^146*x^77 + a^62*x^76 + a^211*x^75 + a^214*x^74 + a^20*x^73 \
							 + a^42*x^72 + a^116*x^71 + a^153*x^70 + a^95*x^69 + a^143*x^68 + a^101*x^67 + a^30*x^66 + a^142*x^65 + a^220*x^64 + a^68*x^63 + a^89*x^62 \
							 + a^111*x^61 + a^180*x^60 + a^109*x^59 + a^204*x^58 + a^159*x^57 + a^87*x^56 + a^54*x^55 + a^193*x^54 + a^187*x^53 + a^193*x^52 + a^162*x^51 \
							 + a^84*x^50 + a^124*x^49 + a^44*x^48 + a*x^47 + a^236*x^46 + a^212*x^45 + a^26*x^44 + a^235*x^43 + a^26*x^42 + a^157*x^41 + a^15*x^40 \
							 + a^203*x^39 + a^223*x^38 + a^132*x^37 + a^191*x^36 + a^103*x^35 + a^118*x^34 + a^52*x^33 + a^27*x^32 + a^81*x^31 + a^64*x^30 + a^11*x^29 \
							 + a^70*x^28 + a^99*x^27 + a^107*x^26 + a^139*x^25 + a^102*x^24 + a^230*x^23 + a^182*x^22 + a*x^21 + a^145*x^20 + a^223*x^19 + a*x^18 \
							 + a^15*x^17 + a^94*x^16 + a^94*x^15 + a^170*x^14 + a^5*x^13 + a^88*x^12 + a^129*x^11 + a^215*x^10 + a^202*x^9 + a^150*x^8 + a^167*x^7 \
							 + a^73*x^6 + a^85*x^5 + a^53*x^4 + a^236*x^3 + a^73*x^2 + a^221*x + a^206")

		self._Sbox_inv = self._P("a^30*x^254 + a^70*x^253 + a^121*x^252 + a^123*x^251 + a^153*x^250 + a^85*x^249 + a^21*x^248 + a^21*x^247 + a^193*x^246 + a^197*x^245 \
							 + a^180*x^244 + a^125*x^243 + a^213*x^242 + a^221*x^241 + a^208*x^240 + a^7*x^239 + a^202*x^238 + a^133*x^237 + a^21*x^236 + a^107*x^235 \
							 + a^9*x^234 + a^239*x^233 + a^136*x^232 + a^29*x^231 + a*x^230 + a^220*x^229 + a^30*x^228 + x^227 + a^142*x^226 + x^225 + a^131*x^224 \
							 + a^144*x^223 + a^199*x^222 + a^124*x^221 + a^178*x^220 + a^51*x^219 + a^27*x^218 + a^249*x^217 + a^10*x^216 + a^227*x^215 + a^124*x^214 \
							 + a^208*x^213 + a^47*x^212 + a^229*x^211 + a^97*x^210 + a^112*x^209 + a^169*x^208 + a^180*x^207 + a^36*x^206 + a^231*x^205 + a^97*x^204 \
							 + a^225*x^203 + a^74*x^202 + a^186*x^201 + a^237*x^200 + a^129*x^199 + a^251*x^198 + a^55*x^197 + a^46*x^196 + a^202*x^195 + a^4*x^194 \
							 + a^171*x^193 + a^112*x^192 + a^26*x^191 + a^201*x^190 + a^189*x^189 + a^35*x^188 + a^237*x^187 + a^168*x^186 + a^238*x^185 + a^111*x^184 \
							 + a^199*x^183 + a^64*x^182 + a^180*x^181 + a^209*x^180 + a^248*x^179 + a^9*x^178 + a^146*x^177 + a^71*x^176 + a^185*x^175 + a^83*x^174 \
							 + a^114*x^173 + a^252*x^172 + a^186*x^171 + a^210*x^170 + a^65*x^169 + a^103*x^168 + a^74*x^167 + a^128*x^166 + a^151*x^165 + a^21*x^164 \
							 + a^104*x^163 + a^79*x^162 + a^135*x^161 + a^134*x^160 + a^116*x^159 + a^32*x^158 + a^115*x^157 + a^15*x^156 + a^173*x^155 + a^82*x^154 \
							 + a^192*x^153 + a^239*x^152 + a^134*x^151 + a^73*x^150 + a^221*x^149 + a^19*x^148 + a^239*x^147 + a^146*x^146 + a^251*x^145 + a^235*x^144 \
							 + a^142*x^143 + a^26*x^142 + a^118*x^141 + a^181*x^140 + a^163*x^139 + a^247*x^138 + a^138*x^137 + a^156*x^136 + a^69*x^135 + a^234*x^134 \
							 + a^253*x^133 + a^220*x^132 + a^28*x^131 + a^143*x^130 + a^6*x^129 + a^132*x^128 + a^66*x^127 + a^40*x^126 + a^126*x^125 + a^74*x^124 \
							 + a^183*x^123 + a^144*x^122 + a^133*x^121 + a^69*x^120 + a^240*x^119 + a^243*x^118 + a^90*x^117 + a^26*x^116 + a^85*x^115 + a^38*x^114 \
							 + a^139*x^113 + a^92*x^112 + a^142*x^111 + a^160*x^110 + a^109*x^109 + a^99*x^108 + a^238*x^107 + a^63*x^106 + a^222*x^105 + a^247*x^104 \
							 + a^242*x^103 + a^85*x^102 + a^21*x^101 + x^100 + a^90*x^99 + a^20*x^98 + a^192*x^97 + a^197*x^96 + a^106*x^95 + a^78*x^94 + a^74*x^93 \
							 + a^85*x^92 + a^218*x^91 + a^141*x^90 + a^144*x^89 + a^164*x^88 + a^198*x^87 + a^134*x^86 + a^18*x^85 + a^155*x^84 + a^231*x^83 + a^59*x^82 \
							 + a^84*x^81 + a^170*x^80 + a^224*x^79 + a^54*x^78 + a^217*x^77 + a^47*x^76 + a^81*x^75 + a^44*x^74 + a^33*x^73 + a^176*x^72 + a^234*x^71 \
							 + a^156*x^70 + a^208*x^69 + a^81*x^68 + a^29*x^67 + a^211*x^66 + a^70*x^65 + a^91*x^64 + a^244*x^63 + a^32*x^62 + a^189*x^61 + a^198*x^60 \
							 + a^104*x^59 + a^127*x^58 + a^164*x^57 + a^118*x^56 + a^229*x^55 + a^43*x^54 + a^40*x^53 + a^65*x^52 + a^114*x^51 + a^252*x^50 + a^84*x^49 \
							 + a^57*x^48 + a^143*x^47 + a^117*x^46 + a^65*x^45 + a^174*x^44 + a^10*x^43 + a^206*x^42 + a^168*x^41 + a^154*x^40 + a^70*x^39 + a^252*x^38 \
							 + a^227*x^37 + a^113*x^36 + a^128*x^35 + a^147*x^34 + a^148*x^33 + a^60*x^32 + a^45*x^31 + a^117*x^30 + a^137*x^29 + a^161*x^28 + a^230*x^27 \
							 + a^136*x^26 + a^86*x^25 + a^67*x^24 + a^200*x^23 + a^178*x^22 + a^34*x^21 + a^250*x^20 + a^254*x^19 + a^92*x^18 + a^170*x^17 + a^10*x^16 \
							 + a^166*x^15 + a^124*x^14 + a^71*x^13 + a^13*x^12 + a^5*x^11 + a^18*x^10 + a^19*x^9 + a^211*x^8 + a^68*x^7 + a^108*x^6 + a^82*x^5 + a^69*x^4 \
							 + a^240*x^3 + a^201*x^2 + a^250*x + a^232")

		# M = matrix(ZZ,16,[
		# 	[0xCF,0x98,0x74,0xBF,0x93,0x8E,0xF2,0xF3,0x0A,0xBF,0xF6,0xA9,0xEA,0x8E,0x4D,0x6E],
		# 	[0x6E,0x20,0xC6,0xDA,0x90,0x48,0x89,0x9C,0xC1,0x64,0xB8,0x2D,0x86,0x44,0xD0,0xA2],
		# 	[0xA2,0xC8,0x87,0x70,0x68,0x43,0x1C,0x2B,0xA1,0x63,0x30,0x6B,0x9F,0x30,0xE3,0x76],
		# 	[0x76,0x33,0x10,0x0C,0x1C,0x11,0xD6,0x6A,0xA6,0xD7,0xF6,0x49,0x07,0x14,0xE8,0x72],
		# 	[0x72,0xF2,0x6B,0xCA,0x20,0xEB,0x02,0xA4,0x8D,0xD4,0xC4,0x01,0x65,0xDD,0x4C,0x6C],
		# 	[0x6C,0x76,0xEC,0x0C,0xC5,0xBC,0xAF,0x6E,0xA3,0xE1,0x90,0x58,0x0E,0x02,0xC3,0x48],
		# 	[0x48,0xD5,0x62,0x17,0x06,0x2D,0xC4,0xE7,0xD5,0xEB,0x99,0x78,0x52,0xF5,0x16,0x7A],
		# 	[0x7A,0xE6,0x4E,0x1A,0xBB,0x2E,0xF1,0xBE,0xD4,0xAF,0x37,0xB1,0xD4,0x2A,0x6E,0xB8],
		# 	[0xB8,0x49,0x87,0x14,0xCB,0x8D,0xAB,0x49,0x09,0x6C,0x2A,0x01,0x60,0x8E,0x4B,0x5D],
		# 	[0x5D,0xD4,0xB8,0x2F,0x8D,0x12,0xEE,0xF6,0x08,0x54,0x0F,0xF3,0x98,0xC8,0x7F,0x27],
		# 	[0x27,0x9F,0xBE,0x68,0x1A,0x7C,0xAD,0xC9,0x84,0x2F,0xEB,0xFE,0xC6,0x48,0xA2,0xBD],
		# 	[0xBD,0x95,0x5E,0x30,0xE9,0x60,0xBF,0x10,0xEF,0x39,0xEC,0x91,0x7F,0x48,0x89,0x10],
		# 	[0x10,0xE9,0xD0,0xD9,0xF3,0x94,0x3D,0xAF,0x7B,0xFF,0x64,0x91,0x52,0xF8,0x0D,0xDD],
		# 	[0xDD,0x99,0x75,0xCA,0x97,0x44,0x5A,0xE0,0x30,0xA6,0x31,0xD3,0xDF,0x48,0x64,0x84],
		# 	[0x84,0x2D,0x74,0x96,0x5D,0x77,0x6F,0xDE,0x54,0xB4,0x8D,0xD1,0x44,0x3C,0xA5,0x94],
		# 	[0x94,0x20,0x85,0x10,0xC2,0xC0,0x01,0xFB,0x01,0xC0,0xC2,0x10,0x85,0x20,0x94,0x01]
		# ])

		# self._G = matrix(self._K,M.nrows())

		# for i in xrange(M.nrows()):
		# 	for j in xrange(M.ncols()):
		# 		self._G[i,j] = self._K.fetch_int(M[i,j])

		if self._key is not None:
			self.__KeyExpansion(self._key)
		else:
			raise ValueError("key must be presented")

	def decrypt(self,ciphertext):
		r'''
			Decryption routine
		'''

		state = self.__to_state(ciphertext)

		state = self.__X(self._rks[9],state)

		for i in xrange(8,-1,-1):
			state = self.__L_inv(state)
			state = self.__S_inv(state)
			state = self.__X(self._rks[i],state)

		return self.__from_state(state)

	def encrypt(self,plaintext):
		r'''
			Encryption routine
		'''

		state = self.__to_state(plaintext)

		state = self.__X(self._rks[0],state)

		for i in xrange(1,10):
			state = self.__S(state)
			state = self.__L(state)
			state = self.__X(self._rks[i],state)

		return self.__from_state(state)

	def __F(self,key,l,r):
		r'''
			The F routine
		'''	

		state = l[:]

		state = self.__X(key,state)
		state = self.__S(state)
		state = self.__L(state)

		state = self.__X(r,state)

		return state,l

	def __from_state(self,state):
		r'''
			Convert the iternal representation to the corresponding integer 
		'''	
		return ZZ([g.integer_representation() for g in state[::-1]],2^8)

	def __KeyExpansion(self,K):
		r'''
			Generate subkeys for all round of encryption based on the given master key
		'''

		C = [self.__L(self.__to_state(i)) for i in xrange(1,33)]

		self._rks[0] = self.__to_state(ZZ(K).digits(2^128,padto=2)[1])
		self._rks[1] = self.__to_state(ZZ(K).digits(2^128,padto=2)[0])

		l = self._rks[0]
		r = self._rks[1]

		for i in xrange(4):
			for j in xrange(8):
				l,r = self.__F(C[i*8+j],l,r)
			self._rks[(i+1)*2]   = l
			self._rks[(i+1)*2+1] = r

	def __l(self,state):
		r'''
			Basic operation l(x)
		'''
		t = self._K(0)

		for i,x in enumerate([148, 32, 133, 16, 194, 192, 1, 251, 1, 192, 194, 16, 133, 32, 148, 1]):
			t += self._K.fetch_int(x)*state[i]

		return t

	def __L(self,state):
		r'''
			The L routine
		'''
		for i in xrange(16):
			state = self.__R(state)

		return state

	def __L_inv(self,state):
		r'''
			The L^{-1} routine
		'''
		for i in xrange(16):
			state = self.__R_inv(state)

		return state

	def __R(self,state):
		r'''
			The R routine
		'''
		return [self.__l(state)] + state[:15]

	def __R_inv(self,state):
		r'''
			The R^{-1} routine
		'''
		return state[1:] + [self.__l(state[1:] + [state[0]])]

	def __S(self,state):
		r'''
			The S routine (very slow)
		'''	
		return [self._Sbox.subs(g) for g in state]

	def __S_inv(self,state):
		r'''
			The S^{-1} routine (very slow)
		'''	
		return [self._Sbox_inv.subs(g) for g in state]

	def selfTesting(self):
		r'''
			Testing the class using data from the specification
		'''

		test_S = [0xffeeddccbbaa99881122334455667700,0xb66cd8887d38e8d77765aeea0c9a7efc,0x559d8dd7bd06cbfe7e7b262523280d39,0x0c3322fed531e4630d80ef5c5a81c50b,0x23ae65633f842d29c5df529c13f5acda]

		# Direct trasformation
		for i in xrange(len(test_S)-1):
			state = self.__to_state(test_S[i])
			state = self.__S(state)

			if (self.__from_state(state) != test_S[i+1]):
				print "{0:032X} ? {1:032X}".format(self.__from_state(state),test_S[i+1])
				raise TypeError("Self testing of '__S' fail")

		# Inverse trasformation
		for i in xrange(len(test_S)-1,0,-1):
			state = self.__to_state(test_S[i])
			state = self.__S_inv(state)

			if (self.__from_state(state) != test_S[i-1]):
				print "{0:032X} ? {1:032X}".format(self.__from_state(state),test_S[i-1])
				raise TypeError("Self testing of '__S_inv' fail")

		test_R = [0x00000000000000000000000000000100,0x94000000000000000000000000000001,0xa5940000000000000000000000000000,0x64a59400000000000000000000000000,0x0d64a594000000000000000000000000]

		# Direct trasformation
		for i in xrange(len(test_R)-1):
			state = self.__to_state(test_R[i])
			state = self.__R(state)

			if (self.__from_state(state) != test_R[i+1]):
				print "{0:032X} ? {1:032X}".format(self.__from_state(state),test_R[i+1])
				raise TypeError("Self testing of '__R' fail")

		# Inverse trasformation
		for i in xrange(len(test_R)-1,0,-1):
			state = self.__to_state(test_R[i])
			state = self.__R_inv(state)

			if (self.__from_state(state) != test_R[i-1]):
				print "{0:032X} ? {1:032X}".format(self.__from_state(state),test_R[i-1])
				raise TypeError("Self testing of '__R_inv' fail")

		test_L = [0x64a59400000000000000000000000000,0xd456584dd0e3e84cc3166e4b7fa2890d,0x79d26221b87b584cd42fbc4ffea5de9a,0x0e93691a0cfc60408b7b68f66b513c13,0xe6a8094fee0aa204fd97bcb0b44b8580]

		# Direct trasformation
		for i in xrange(len(test_L)-1):
			state = self.__to_state(test_L[i])
			state = self.__L(state)

			if (self.__from_state(state) != test_L[i+1]):
				print "{0:032X} ? {1:032X}".format(self.__from_state(state),test_L[i+1])
				raise TypeError("Self testing of '__L' fail")

		# Inverse trasformation
		for i in xrange(len(test_L)-1,0,-1):
			state = self.__to_state(test_L[i])
			state = self.__L_inv(state)

			if (self.__from_state(state) != test_L[i-1]):
				print "{0:032X} ? {1:032X}".format(self.__from_state(state),test_L[i-1])
				raise TypeError("Self testing of '__L_inv' fail")

		test_K = 0x8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef

		test_rks = [0x8899aabbccddeeff0011223344556677, 0xfedcba98765432100123456789abcdef, 0xdb31485315694343228d6aef8cc78c44, 0x3d4553d8e9cfec6815ebadc40a9ffd04, 
			 		0x57646468c44a5e28d3e59246f429f1ac, 0xbd079435165c6432b532e82834da581b, 0x51e640757e8745de705727265a0098b1, 0x5a7925017b9fdd3ed72a91a22286f984, 
			 		0xbb44e25378c73123a5f32f73cdb6e517, 0x72e9dd7416bcf45b755dbaa88e4a4043]

		rks = self._rks[:]
		self.__KeyExpansion(test_K)

		if ([self.__from_state(k) for k in self._rks] != test_rks):
			raise TypeError("Self testing of KeyExpansion fail")

		M = 0x1122334455667700ffeeddccbbaa9988
		C = 0x7f679d90bebc24305a468d42b9d4edcd

		if (self.encrypt(M) != C):
			raise TypeError("Self testing of encryption fail")

		if (self.decrypt(C) != M):
			raise TypeError("Self testing of decryption fail")		

		self._rks = rks[:]

	def __to_state(self,state):
		r'''
			Convert an integer to the iternal representation
		'''	
		return [self._K.fetch_int(g) for g in ZZ(state).digits(base=2^8,padto=16)[::-1]]

	def __X(self,key,state):
		r'''
			The X routine
		'''
		return [g[0]+g[1] for g in zip(key,state)]
