import sys
import genedist
from os.path import exists,isfile,isdir

PROJECT_DIR = "/Users/lizkennedy/Dropbox/CurrentProj/Thesis/"

#checks if a file exists
#should_exist indicates whether we expect the file to exist or not
#if it doesnt match our expectation, it will stop the program
def check_file(path,should_exist=True):
	if not should_exist and exists(path):
		print "Cant Create File, Already Exists: %s" % path
		sys.exit(-1)
	elif should_exist and (not exists(path) or not isfile(path)):
		print "File Not Found: %s" % path
		sys.exit(-1)

#checks that the directory exists
def check_dir(path):
	if not exists(path) or not isdir(path):
		print "Directory Not Found: %s" % path

CPG_FILENAME = PROJECT_DIR + "MESA_CpGLoc_4_pylmm.txt"
print "Cpg locations filename"
print "CPG_FILENAME = %s\n" % (CPG_FILENAME)
check_file(CPG_FILENAME)

GENE_REF_FILE = PROJECT_DIR + "MESA_geneLoc_4_pylmm.txt"
print "Gene reference file"
print "GENE_REF_FILE = %s\n" % (GENE_REF_FILE)
check_file(GENE_REF_FILE)
genedist.load_file(GENE_REF_FILE) #load the genes into the database

CLOSEST_GENE_FILE = PROJECT_DIR + "MESA_Closest_Gene.txt"
print "Closest gene file"
print "CLOSEST_GENE_FILE = %s\n" % (CLOSEST_GENE_FILE)
check_file(CLOSEST_GENE_FILE,False)

with open(CPG_FILENAME) as infile:
	with open(CLOSEST_GENE_FILE,'w') as outfile:
		outfile.write("CPG_ID\tABS_DIST\n")
		for line in infile:
			cpgnm,chrm,pos = line.strip().split("\t")
			abs_dist = genedist.find_abs_closest(cpgnm,chrm,pos)
			print "CpG: " + str(chrm)
			if abs_dist is None:
				print "ERROR: No genes found for %s" % (cpgnm)
			else:
				outfile.write("%s\t%s\n" % (cpgnm,abs_dist))
