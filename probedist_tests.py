import unittest
import probedist
import os

class TestProbeDist(unittest.TestCase):
	def test_format_category(self):
		self.assertEquals(probedist.format_category("PREF"),0)
		self.assertEquals(probedist.format_category("PREFERRED"),0)
		self.assertEquals(probedist.format_category("pref"),0)
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

	def test_within_gene(self):
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("LARGE","CHR1",2999)
		self.assertEquals(priority,0)
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST","CHR3",35)
		self.assertEquals(priority,3)
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST","CHR2",100)
		self.assertEquals(priority,0)
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST","CHR2",200)
		self.assertEquals(priority,0)

	def test_probedist(self):
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST","CHR1",100)
		self.assertEquals(start,45)
		self.assertEquals(stop,10)
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST","CHR2",10)
		self.assertEquals(start,100)
		self.assertEquals(stop,200)

	def test_prioritize_ingene(self):
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST1","CHR1",24)
		self.assertEquals(start,32)
		self.assertEquals(stop,21)
		self.assertEquals(priority,0)
		self.assertEquals(category,2)

	def test_prioritize_ingene_multi(self):
		#probes should be prioritized by category first
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST2","CHR10",320)
		self.assertEquals(start,200)
		self.assertEquals(stop,400)
		self.assertEquals(category,0)
		#then by start distance
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST2","CHR10",490)
		self.assertEquals(start,540)
		self.assertEquals(stop,320)
		self.assertEquals(category,1)

	def test_higher_priority_lower_category(self):
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST3","CHR3",175)
		self.assertEquals(start,372)
		self.assertEquals(stop,125)
		self.assertEquals(strand,'-')
		self.assertEquals(category,1)

	def test_start_dist(self):
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST1","CHR1",24)
		self.assertEquals(start_dist,8)
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST1","CHR1",98)
		self.assertEquals(start_dist,2)
		prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST1","CHR2",98)
		self.assertEquals(start_dist,None)

	def test_random_diff_chrm(self):
		results = set([])
		#someone needs to tell me the probability of randint performing a coin flip and landing on heads 100 times in a row
		for i in range(0,100):
			prbchrm,start,stop,strand,start_dist,category,priority = probedist.closest_probe("TEST","CHR3",98)
			results.add(start)
		self.assertTrue(len(results),1)
		self.assertTrue(45 in results)

def loadtestdata():
	TEST_FILENAME = "test.txt"
	with open(TEST_FILENAME,"w") as outfile:
		outfile.write("TEST\tCHR1\t10\t45\t-\tPREF\n")
		outfile.write("LARGE\tCHR1\t1\t3000\t+\tOTH\n")
		outfile.write("TEST\tCHR2\t100\t200\t+\tOTH\n")
		outfile.write("TEST1\tCHR1\t21\t32\t-\tOTH\n")
		outfile.write("TEST1\tCHR1\t100\t200\t+\tPREF\n")
		outfile.write("TEST2\tCHR10\t200\t400\t+\tPREF\n")
		outfile.write("TEST2\tCHR10\t300\t500\t+\tSEC\n")
		outfile.write("TEST2\tCHR10\t320\t540\t-\tSEC\n")
		outfile.write("TEST3\tCHR3\t125\t372\t-\tSEC\n")
		outfile.write("TEST3\tCHR3\t200\t400\t+\tPREF\n")
	probedist.loadprobelocs(TEST_FILENAME)
	os.remove(TEST_FILENAME)

if __name__ == '__main__':
	loadtestdata()
	unittest.main()
