import unittest
import genedist
import os

class TestGeneDist(unittest.TestCase):
	def test_count(self):
		genes = genedist.stmt("SELECT * FROM genes")
		self.assertEquals(len(genes),1)

	def test_findbetween(self):
		genes = genedist.findbetween("CHR1",2,27)
		self.assertEquals(len(genes),1)


def loadtestdata():
	TEST_FILENAME = "test.txt"
	with open(TEST_FILENAME,"w") as outfile:
		outfile.write("GENE\tCHRM\tSTRAND\tSTART\tSTOP\n")
		outfile.write("TEST1\tCHR1\t+\t20\t30\n")
	genedist.load_file(TEST_FILENAME)
	os.remove(TEST_FILENAME)

if __name__ == '__main__':
	loadtestdata()
	unittest.main()
