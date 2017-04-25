import unittest
import probedist
import os

class TestProbeDist(unittest.TestCase):
	def test_format_category(self):
		self.assertEquals(probedist.format_category("PRIM"),0)
		self.assertEquals(probedist.format_category("PRIMARY"),0)
		self.assertEquals(probedist.format_category("pri"),0)
		self.assertEquals(probedist.format_category("sec"),1)
		self.assertEquals(probedist.format_category("second"),1)
		self.assertEquals(probedist.format_category("secondary"),1)
		self.assertEquals(probedist.format_category("oth"),2)
		self.assertEquals(probedist.format_category("OTHER"),2)
	def test_frmt_chrm(self):
		self.assertEquals(probedist.frmt_chrm("chr15"),"15")
		self.assertEquals(probedist.frmt_chrm("15"),"15")
		self.assertEquals(probedist.frmt_chrm("CHR1"),"1")
		self.assertEquals(probedist.frmt_chrm("1"),"1")


	def test_probedist(self):
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST","CHR1",100)
		self.assertEquals(start,45)
		self.assertEquals(stop,10)
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST","CHR2",10)
		self.assertEquals(start,100)
		self.assertEquals(stop,200)

	def test_prioritize_ingene(self):
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST1","CHR1",24)
		self.assertEquals(start,21)
		self.assertEquals(priority,0)
		self.assertEquals(category,2)

	def test_start_dist(self):
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST1","CHR1",24)
		self.assertEquals(start_dist,3)
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST1","CHR1",98)
		self.assertEquals(start_dist,2)

			

if __name__ == '__main__':
	TEST_FILENAME = "test.txt"
	with open(TEST_FILENAME,"w") as outfile:
		outfile.write("TEST\tCHR1\t10\t45\t-\tPRI\n")
		outfile.write("TEST\tCHR2\t100\t200\t+\tOTH\n")
		outfile.write("TEST1\tCHR1\t32\t21\t-\tOTH\n")
		outfile.write("TEST1\tCHR1\t100\t200\t+\tPRI\n")
	probedist.loadprobelocs(TEST_FILENAME)
	unittest.main()
	os.remove(TEST_FILENAME)
