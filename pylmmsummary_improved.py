import sys
import gzip
import genedist
from os.path import exists,isfile,isdir

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
def load_skipgenes(fname):
	genes=set()
	with open(fname) as infile:
		for line in infile:
			genes.add(line.strip())
	return genes
#DROPLINES_PATH = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/infiles/tossIDXs/"
DROPLINES_PATH = "/home/ekennedy/Projects/Thesis/input/tossIDXs/"
print "Directory that contains row indicies to drop"
print "DROPLINES_PATH = %s\n" % (DROPLINES_PATH)
check_dir(DROPLINES_PATH)

START_ID = 18001
STOP_ID = 19445


#PROBE_LOCATIONS_FILENAME = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/infiles/MESA_geneLoc_4_pylmm.txt"
PROBE_LOCATIONS_FILENAME = "/home/ekennedy/Projects/Thesis/input/GTP_geneLoc_4_pylmm.txt"
print "Probe locations filename"
print "PROBE_LOCATIONS_FILENAME = %s\n" % (PROBE_LOCATIONS_FILENAME)
check_file(PROBE_LOCATIONS_FILENAME)

#CPG_FILENAME = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/infiles/MESA_CpGLoc_4_pylmm.txt"
CPG_FILENAME = "/home/ekennedy/Projects/Thesis/input/GTP_CpGLoc_4_pylmm.txt"
print "Cpg locations filename"
print "CPG_FILENAME = %s\n" % (CPG_FILENAME)
check_file(CPG_FILENAME)

#ROWNAME_PATH = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/infiles/MESA_exp_rownames.txt"
ROWNAME_PATH = "/home/ekennedy/Projects/Thesis/input/"
print "'File or directory containing rownames for genes(?)"
print "ROWNAMES_PATH = %s\n" % (ROWNAME_PATH)

#PYLMM_PATH = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/outfiles/newnames/"
PYLMM_PATH = "/home/ekennedy/Projects/Thesis/output/"
print "Directory containing outputs from pylmm"
print "PYLMM_PATH = %s\n" % (PYLMM_PATH)
check_dir(PYLMM_PATH)

#CPG_COUNTS_FILE = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/outfiles/MESA_pylmm_summary_report_nCGI_neg11_"+ str(STOP_ID) + ".csv"
CPG_COUNTS_FILE = "/home/ekennedy/Projects/Thesis/output/GTP_pylmm_summary_report" + ".csv"
print "Pylmm CPG summary file"
print "CPG_COUNTS_FILE = %s\n" % (CPG_COUNTS_FILE)
check_file(CPG_COUNTS_FILE + str(STOP_ID) + ".csv",False)

#CLOSER_GENES_FILE = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/outfiles/MESA_closer_nCGI_neg11_" + str(STOP_ID) + ".csv"
CLOSER_GENES_FILE = "/home/ekennedy/Projects/Thesis/output/GTP_closer" + ".csv"
print "File containing genes which are closer to the cpg than the pylmm significant gene"
print "CLOSER_GENES_FILE = %s\n" % (CLOSER_GENES_FILE)
check_file(CLOSER_GENES_FILE + str(STOP_ID) + ".csv",False)

#GENE_REF_FILE = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/infiles/AVGrefSeq.txt"
GENE_REF_FILE = "/home/ekennedy/Projects/Thesis/input/AVGrefSeq.txt"
print "Gene reference file"
print "GENE_REF_FILE = %s\n" % (GENE_REF_FILE)
check_file(GENE_REF_FILE)
genedist.load_file(GENE_REF_FILE) #load the genes into the database

THRESHOLD = pow(10,-11)
print "THRESHOLD = %d\n" % (THRESHOLD)

WEAK_THRESHOLD = pow(10,-5)
print "WEAK THRESHOLD = %d\n" % (WEAK_THRESHOLD)

#SKIP_GENES = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/infiles/MESA_skip_genes_nCGI.txt"
SKIP_GENES = None
if SKIP_GENES is not None:
	print "File containing expression probe file names to skip"
	print "SKIP_GENES = %s\n" % (SKIP_GENES)
	check_file(SKIP_GENES)
	skip_genes=load_skipgenes(SKIP_GENES)
else:
	skip_genes=set()

#set the appropriate naming convetion to determine how to iterate
#through the files
NAMING_CONVENTION = "GTP"
if NAMING_CONVENTION not in ["GTP","MESA"]:
	print "Invalid naming convention. Please use either MESA, or GTP."
	sys.exit(-1)
else:
	print "NAMING_CONVENTION = %s\n" % (NAMING_CONVENTION)

if NAMING_CONVENTION == "MESA":
	check_file(ROWNAME_PATH)	
else:
	check_dir(ROWNAME_PATH)
def gtp_rowname(prefix,idx):
	return ROWNAME_PATH + "%s%d_rownames.txt" % (prefix,idx)


#load all the probes from the correct rownames file
def loadprobenames(fname):
	probes = []
	with open(fname) as infile:
		for line in infile:
			for prb in line.strip().split(","):
				probes.append(prb)
	return probes

def gtp_name(prefix,row,col):
	return PYLMM_PATH + "%s%d_%d" % (prefix,row,col)

def mesa_name(idx):
	return PYLMM_PATH + "%d" % (idx)

def loaddroplines(lineno):
	lines = set()
	try:
		with open(DROPLINES_PATH + "probe%s.csv" % (lineno)) as infile:
			for line in infile:
				lines.add(int(float(line.strip())))
	except IOError:
		pass
	return lines

#call the appropriate function to calculate the distance between the position and the cpg strand
#based on the direction of the strand (+ or -)
def calcdist(pos,start,stop,strand):
        if strand == "+":
                return calcposdist(pos,start,stop)
        else:
                return calcnegdist(pos,start,stop)

def calcposdist(pos,start,stop):
        if pos < start:
                return pos - start
        elif pos > stop:
                return pos - stop
        else:
                return 0

def calcnegdist(pos,start,stop):
        if pos > start:
                return start - pos
        elif pos < stop:
                return stop - pos
        else:
                return 0

def intbucket(val):
	if val == 0:
		return 0
	elif val >= 2000:
		return val / 1000 + 19
	else:
		return val / 100 + 1

def bucket(val):
	if val >= 0: return intbucket(val)
	else: return -1 * intbucket(-val)

def loadprobelocs():
        probes = {}
        with open(PROBE_LOCATIONS_FILENAME) as infile:
                lineno = 0
                for line in infile:
                        lineno += 1
                        nm,chrm,start,stop,strand = line.strip().split("\t")
                        #if the strand is negative we need to swap the start and stop locations
                        if strand == "-":
                                tmp = stop
                                stop = start
                                start = tmp
                        if start != "NA" and stop != "NA":
				probes[nm] = (chrm,int(start),int(stop),strand,lineno)
        return probes

def loadcpglocs():
	cpglocs = []
	with open(CPG_FILENAME) as infile:
		for line in infile:
			cpgnm,chrm,pos = line.strip().split("\t")
			#cpg,chrm,position,bonferroni significant,distal,trans,weak,buckets
			cpglocs.append([cpgnm,chrm,pos,0,0,0,0,{}])
	return cpglocs 

#file_handle = the file object to write entries to
#sig_gene = the significant gene that the genes in the list are closer to the cpg
#cpgname = the cpg that the genes are close to
#cpgpos = the position of the cpg position
#loc_type = How close the genes in the gene list are relative to the significant gene to the cpg. BETWEEN,TSS, or CLOSER
#gene_list = list of gnees that are closer to the cpg than the significant gene
def write_closer_genes(file_handle,sig_gene,cpgname,cpgpos,dist,loc_type,gene_list):
	for gene,start,stop,strand in gene_list:
		writetofile(file_handle,[sig_gene,cpgname,dist,loc_type,gene,calcdist(cpgpos,start,stop,strand)])

def writetofile(file_handle,entries):
	file_handle.write("%s\n" % (",".join([str(x) for x in entries])))

def get_fname(path):
	return  path.split("/")[-1]

def open_file(fname):
	if NAMING_CONVENTION == "MESA":
		return gzip.open(fname + ".gz")
	else:
		return open(fname)

def loadpylmmresults(gene_fname,probename):
	if probename not in probelocs.keys():
		print("Probename not in ProbeLOCS: " + probename)
		return
	if get_fname(gene_fname) in skip_genes:
		return
	prbchrm,start,stop,strand,prbline = probelocs[probename]
	droplines = loaddroplines(prbline)
	#SNP_ID	BETA	BETA_SD	F_STAT	P_VALUE
	with open_file(gene_fname) as infile:
		infile.readline() # drop the header
		lineno = 0
		for line in infile: #iterate through gene file line by line
			lineno += 1

			#if we've already determined the relationship of this gene 
			#to this cpg isn't significant skip it
			if lineno in droplines:
				continue	
			
			#load the gene
			snp,beta,beta_sd,fstat,pval = line.strip().split("\t")

			#find the cpg for this particular gene (?) I'm not sure why we're using the same
			#indexing as the genes. I imagine there is a cpg listing for each gene, not positive though
			cpgname,chrm,pos,bs,distal,trans,weak,buckets = cpglocations[lineno-1]
			
			if pos != "NA" and pval != "NA":
				if float(pval) < WEAK_THRESHOLD:# and float(fstat)>0:
					if float(pval) < THRESHOLD:
						pos = int(pos)
						bs += 1
						if chrm != prbchrm:
							trans += 1
							writetofile(closer_file,[probename,cpgname,"NA","TRANS","NA","NA"])
						else:
							dist = calcdist(pos,start,stop,strand)
							if dist == 0:
								writetofile(closer_file,[probename,cpgname,0,"IN","NA","NA"])
								bucket_id = bucket(dist)
								buckets[bucket_id] = buckets.get(bucket_id,0) + 1

							else:
								if abs(dist) > 50000:
									distal += 1
									writetofile(closer_file,[probename,cpgname,dist,"DISTAL","NA","NA"])
								else:

									#find any genes with a closer tss
									tss = genedist.findtss(chrm,start,pos)
									#find any genes with a closer absolute distance to the cpg
									closer = genedist.findcloser(chrm,dist,pos)
									#find any genes between the significant gene and cpg
									between = genedist.findbetween(chrm,start,pos)

									bucket_id = bucket(dist)
									buckets[bucket_id] = buckets.get(bucket_id,0) + 1
									#format is sig gene, cpg,closer gene, type,distance
									#if this is the closest gene to the cpg, it means all the other "closer" buckets will be edmpty
									if len(tss) == len(closer) == len(between) == 0:
										writetofile(closer_file,[probename,cpgname,dist,"CLOSEST","NA","NA"])
									else:
										writefn = lambda pos_type,genes: write_closer_genes(closer_file,probename,cpgname,pos,dist,pos_type,genes)
										writefn("BETWEEN",between)
										writefn("TSS",tss)
										writefn("CLOSER",closer)
					else:
						weak += 1 #if it met the weak threshold but not the strong, increment the weak threshold counter
					#save the updated counts for the cpg back into the cpg locations array
					cpglocations[lineno-1] = [cpgname,chrm,pos,bs,distal,trans,weak,buckets]

#naming conventions for different runs of pylmm
#all files are named in the format of prefix_group_gene for input files
#and pylmm_group_gene for output files from pylmm
#where prefix, group, gene, and pylmm are all specified for gene groups below
probe_group_names = [{
	"prefix":"HGCC_test_y",
	"groups":[1],
	"genes":range(10),
	"pylmm":"HGCCtest_"
	},
	{
	"prefix":"HGCC_test_y",
	"groups":range(3,9),
	"genes":range(20),
	"pylmm":"HGCCtest_"
	},
	{
	"prefix":"pyLMM_y",
	"groups":range(9,276),
	"genes":range(0,50),
	"pylmm":"pyLMM_y"
	},
	{
	"prefix":"pyLMM_y",
	"groups":[276],
	"genes":range(0,23),
	"pylmm":"pyLMM_y"
	}]

#loads the probe locations
print "Loading probe locations"
probelocs = loadprobelocs()

#loads the cpg locations
print "Loading cpg locations"
cpglocations = loadcpglocs()

#MAIN LOGIC HERE
#for each of the probe groups (naming conventions for differnt sets of data)
#read through and find significant associations between genes and cpgs, and count
#distance information for each cpg
#also record any closer genes to the CPG than the significant gene into closer_file
closer_file = open(CLOSER_GENES_FILE,"w")
if NAMING_CONVENTION == "GTP":
	for grp in probe_group_names:
		for i in grp["groups"]:
			print "Loading gene names for %s %d" % (grp["prefix"],i)
			probenames = loadprobenames(gtp_rowname(grp["prefix"],i))
			for j in grp["genes"]:
				print "running " + str(grp["prefix"]) + " " + str(i) + " " + str(j)
				loadpylmmresults(gtp_name(grp["pylmm"],i,j),probenames[j])
elif NAMING_CONVENTION == "MESA":
	probenames = loadprobenames(ROWNAME_PATH)
	for idx in range(START_ID,STOP_ID+1):
		print "running " + str(mesa_name(idx+1)) + ": " + str(probenames[idx])
		loadpylmmresults(mesa_name(idx+1),probenames[idx])	
closer_file.close()

#WRITE THE CPG DISTANCE COUNTS OUT TO A FILE
with open(CPG_COUNTS_FILE,"w") as outfile:
	buckets = set()
	for entry in cpglocations:
		buckets = buckets.union(entry[-1].keys())
	buckets = sorted(buckets)
	outfile.write("cpgname,chrm,position,bs,distal,trans,weak," + ",".join([str(x) for x in buckets]) + "\n")
	for entry in cpglocations:
		cpgname,chrm,pos,bs,distal,trans,weak,cnts = entry
		bucket_cnts = ",".join([str(cnts.get(x,0)) for x in buckets])
		outfile.write("%s,%s,%s,%d,%d,%d,%d,%s\n" % (cpgname,chrm,pos,bs,distal,trans,weak,bucket_cnts))
