
import sys
import gzip
import genedist
import probedist
from os.path import exists,isfile,isdir

#when the probe locs file is loaded, there is a categorical column indiciating if its the preferred, secondary, or other  probe location type.
#there will be multiple locations for each probe included in probe locs

#when making comparisons between distances between cpg & probe locs apply the following logic
#if a cpg is within the corresponding probe loc, choose it
  #if multiple, order by primary / secondary / other
    #if multiple in the same category, if the start sites are within 500bp of each other and one has a gene name not starting in "ILMN" - otherwise-choose closest to start site

#if a cpg is within 1500 base pairs upstream of the start site, choose it
  #if multiple, order by primary / secondary / other
    #if multiple in the same category, choose closest to start site

#if a cpg is on the same chrm
  #if multiple, order by primary / secondary / other
    #if multiple in the same category, choose closest to start site

#if a cpg is on the same chrm
  #if multiple, order by primary / secondary / other
    #if multiple in the same category, choose closest to start site

#if a cpg is om a different chrm
  #if multiple, order by primary / secondary / other
    #if multiple, choose the first one

#include the probe information in the output.  probe id, chrm, start, end

#include distance to start site from cpg

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
### WE NEED TO GET RID OF DROPLINES ALTOGETHER ###
#DROPLINES_PATH = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/infiles/tossIDXs/"
#DROPLINES_PATH = "/home/ekennedy/Projects/Thesis/input/tossIDXs/
#print "Directory that contains row indicies to drop"
#print "DROPLINES_PATH = %s\n" % (DROPLINES_PATH)
#check_dir(DROPLINES_PATH)

START_ID = 18001
STOP_ID = 19445


PROBE_LOCATIONS_FILENAME = "/Users/lizkennedy/Dropbox/CurrentProj/Thesis/exp_geneLoc_4_pylmm.txt"
print "Probe locations filename"
print "PROBE_LOCATIONS_FILENAME = %s\n" % (PROBE_LOCATIONS_FILENAME)
check_file(PROBE_LOCATIONS_FILENAME)

CPG_FILENAME = "/Users/lizkennedy/Dropbox/CurrentProj/Thesis/MESA_CpGLoc_4_pylmm.txt"
#CPG_FILENAME = "/Users/lizkennedy/Dropbox/CurrentProj/Thesis/GTP_CpGLoc_4_pylmm.txt"
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


CLOSER_INPUT_FILE = "/Users/lizkennedy/Dropbox/CurrentProj/Thesis/MESA_closer2.txt"

#CLOSER_GENES_FILE = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/outfiles/MESA_closer_nCGI_neg11_" + str(STOP_ID) + ".csv"
CLOSER_GENES_FILE = "/Users/lizkennedy/Dropbox/CurrentProj/Thesis/MESA_adj_closer.csv"
print "File containing genes which are closer to the cpg than the pylmm significant gene"
print "CLOSER_GENES_FILE = %s\n" % (CLOSER_GENES_FILE)
check_file(CLOSER_GENES_FILE + str(STOP_ID) + ".csv",False)

#GENE_REF_FILE = "/gpfs/pace1/home/het1/emkenn2-emory/data/MESA/infiles/AVGrefSeq.txt"
GENE_REF_FILE = "/Users/lizkennedy/Dropbox/CurrentProj/Thesis/AVGrefSeq.txt"
print "Gene reference file"
print "GENE_REF_FILE = %s\n" % (GENE_REF_FILE)
check_file(GENE_REF_FILE)
genedist.load_file(GENE_REF_FILE) #load the genes into the databse

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
#Liz made changes here to try to fix a typeERROR*****************
def calcdist(pos,start,stop,strand):
	if strand == "+":
		return calcposdist(int(pos),int(start),int(stop))
        else:
		return calcnegdist(int(pos),int(start),int(stop))

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

#THIS FUNCTION HAS BEEN MIGRATED TO PROBEDIST AND IS NO LONGER USED
#def loadprobelocs():
#        probes = {}
#        with open(probe_locations_filename) as infile:
#                lineno = 0
#                for line in infile:
#                        lineno += 1
#                        nm,chrm,start,stop,strand = line.strip().split("\t")
#                        #if the strand is negative we need to swap the start and stop locations
#                       if strand == "-":
#                                tmp = stop
#                                stop = start
#                                start = tmp
#                        if start != "NA" and stop != "NA":
#				probes[nm] = (chrm,int(start),int(stop),strand,lineno)
#        return probes

def loadcpglocs():
	cpglocs = {}
	with open(CPG_FILENAME) as infile:
		for line in infile:
			cpgnm,chrm,pos = line.strip().split("\t")
			cpglocs[cpgnm]=[chrm,pos]
	return cpglocs 

#file_handle = the file object to write entries to
#sig_gene = the significant gene that the genes in the list are closer to the cpg
#cpgname = the cpg that the genes are close to
#cpgpos = the position of the cpg position
#loc_type = How close the genes in the gene list are relative to the significant gene to the cpg. BETWEEN,TSS, or CLOSER
#gene_list = list of gnees that are closer to the cpg than the significant gen
def write_closer_genes(file_handle,cpg,cpgpos,probe,prbchrm,prb_start,prb_stop,prb_strand,start_dist,sig_gene,sig_dist,loc_type,gene_list,pval,fstat,beta,beta_sd):
	for gene,start,stop,strand in gene_list:
		oth_dist = calcdist(cpgpos,start,stop,strand)
		writetofile(file_handle,[cpg,probe,prbchrm,prb_start,prb_stop,prb_strand,start_dist,sig_gene,sig_dist,loc_type,gene,oth_dist,pval,fstat,beta,beta_sd])

def writetofile(file_handle,entries):
	file_handle.write("%s\n" % ("\t".join([str(x) for x in entries])))

def loadcloser(filename):
	with open(filename) as infile:
		infile.readline()#drop header line
		for line in infile:
			yield line.strip().split("\t")

def findcloser(chrm,start,dist,cpgpos):
	tss = genedist.findtss(chrm,start,cpgpos)
	closer = genedist.findcloser(chrm,dist,cpgpos)
	between = genedist.findbetween(chrm,start,cpgpos)
	tss = addstatus("TSS",tss)
	closer = addstatus("CLOSER",closer)
	between = addstatus("BETWEEN",between)
	return between + tss + closer

def addstatus(status,genes):
	coll = []
	for gene,start,stop,strand in genes:
		coll.append((gene,start,stop,strand,status))
	return coll

def findclosest():
	for entry in loadcloser(CLOSER_INPUT_FILE):
		cpg,probe,junk,distance,status,other_gene,other_distance,pval,fstat,beta,beta_sd = entry
		#output format should be
		#cpg,probe,probe_chrm,probe_strand,probe_start,probe_end,start_dist,gene,distance,status,other_gene,other_distance,pval,fstat,beta,beta_sd
		chrm,cpgpos = cpglocations[cpg]
		probeinfo = probedist.closest_probe(probe,chrm,cpgpos)
		if not probeinfo:
			writetofile(closer_file,[cpg,probe,"NA","NA","NA","NA","NA","NA","NA","NOPROBELOC","NA","NA",pval,fstat,beta,beta_sd])
		else:
			prbchrm,prb_start,prb_end,strand,gene,start,stop,start_dist,category,priority = probeinfo
			if not start_dist:
				start_dist = "NA"
			if chrm != prbchrm:
				writetofile(closer_file,[cpg,probe,prbchrm,prb_start,prb_end,strand,start_dist,gene,"NA","TRANS","NA","NA",pval,fstat,beta,beta_sd])
			else:
				dist = calcdist(cpgpos,start,stop,strand)
				#if the distance indicates the cpg falls within the probe, or significantly far away dont bother looking for other closer genes
				if dist == 0:
					writetofile(closer_file,[cpg,probe,prbchrm,prb_start,prb_end,strand,start_dist,gene,0,"IN","NA","NA",pval,fstat,beta,beta_sd])
				elif abs(dist) > 50000:
					writetofile(closer_file,[cpg,probe,prbchrm,prb_start,prb_end,strand,start_dist,gene,dist,"DISTAL","NA","NA",pval,fstat,beta,beta_sd])
				else:
					closer = findcloser(chrm,start,dist,cpgpos)
					#if this is the closest gene to the cpg, it means all the other "closer" buckets will be edmpty
					if len(closer) == 0:
						writetofile(closer_file,[cpg,probe,prbchrm,prb_start,prb_end,strand,start_dist,gene,dist,"CLOSEST","NA","NA",pval,fstat,beta,beta_sd])
					else:
						for oth_gene,oth_start,oth_stop,oth_strand,oth_status in closer:
							oth_dist = calcdist(cpgpos,oth_start,oth_stop,oth_strand)
							writetofile(closer_file,[cpg,probe,prbchrm,prb_start,prb_end,strand,start_dist,gene,dist,oth_status,oth_gene,oth_dist,pval,fstat,beta,beta_sd])	

#naming conventions for different runs of pylmm
#all files are named in the format of prefix_group_gene for input files
#and pylmm_group_gene for output files from pylmm
#where prefix, group, gene, and pylmm are all specified for gene groups below
#loads the probe locations
print "Loading probe locations"
probedist.loadprobelocs(PROBE_LOCATIONS_FILENAME)
#probelocs = loadprobelocs()

#loads the cpg locations
print "Loading cpg locations"
cpglocations = loadcpglocs()

#MAIN LOGIC HERE
#for each of the probe groups (naming conventions for differnt sets of data)
#read through and find significant associations between genes and cpgs, and count
#distance information for each cpg
#also record any closer genes to the CPG than the significant gene into closer_file
closer_file =open(CLOSER_GENES_FILE,"w")
closer_file.write("cpg\tprobe\tprobe_chrm\tprobe_start\tprobe_stop\tprobe_strand\tstart_dist\tgene\tdistance\tstatus\tother_gene\tother_distance\tpval\tfstat\tbeta\tbeta_sd\n")
findclosest()
closer_file.close()
