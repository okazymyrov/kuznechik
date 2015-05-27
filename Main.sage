#!/usr/bin/env sage

r"""
    An example of how to use the class "Kuznechik" from the external code

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

load ./Kuznechik_2_8.sage
load ./Kuznechik.sage
load ./Kuznechik_table.sage
load ./Tools.sage

def main(argv=None):
	K = ZZ.random_element(2^256)
	M = ZZ.random_element(2^128)

	# generate_tables(filename="Tables.sage")
	# return

	# test_performance()
	# return

	# cipher = Kuznechik_2_8(key=K)
	# cipher = Kuznechik(key=K)
	cipher = Kuznechik_table(key=K)

	t1 = cputime()

	cipher.selfTesting()

	C = cipher.encrypt(M)

	if cipher.decrypt(C) == M:
		print "Encryption/Decryption: Pass"
	else:
		print "Encryption/Decryption: Fail"

	t2 = cputime()

	print "====="
	print "Time = {0}".format(t2-t1)


if __name__ == "__main__":
    sys.exit(main())
