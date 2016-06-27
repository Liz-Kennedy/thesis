from math import log
import sqlite3
import json
import subprocess
import sys

CLOSER_FILE = "MESA_closer_neg11.txt"  #holds probe, cpg, distance (and more)

CPG_POSITION_FILE = "MESA_CpGLoc_4_pylmm.txt"

GENE_FILE = "./TARDIS_UPLOAD/AVGrefSeq.txt"

PROBE_TO_GENE_MAP = "Illumina_HT12_genenames.txt"

NUM_BUCKETS = 10

GENE_LENGTHS_FILE = "MESA_circos_links.txt"

def loadprbgenemap(fname):
	amap = {}
	with open(fname) as infile:
		lines = 0
		for line in infile:
			lines += 1
			try:
				prb,gene = line.upper().strip().split(" ")
				amap[prb] = gene
			except:
				print "Failed to parse %d" % (lines)
	return amap
probemap = loadprbgenemap(PROBE_TO_GENE_MAP)

def loadgenes(fname):
	genes = {}
	with open(fname) as infile:
		infile.readline()
		for line in infile:
			nm,chrm,strand,start,end,refseqids = line.upper().strip().split("\t")
			genes[nm] = (int(start),int(end),strand,chrm)
	return genes

genenames = loadgenes(GENE_FILE)

def loadcpgs(fname):
	cpgs = {}
	with open(fname) as infile:
		for line in infile:
			cpgname,chrm,cpgpos = line.upper().strip().split("\t")
			if cpgpos != "NA":
				cpgs[cpgname] = (chrm,int(cpgpos))
	return cpgs

cpg_positions = loadcpgs(CPG_POSITION_FILE)

missingcpg = set()
missingprobe = set()
missinggene = set()

def is_missing(key,amap,missing_entries,warning):
	if key not in amap:
		if key not in missing_entries:
			missing_entries.add(key)
			print warning + " " + key
		return True
	else:
		return False

def formatchrome(chrm):
	return "hs" + chrm[3:]


#keep a set of all the unique cpg, gene pairs
sig_pairs = {}

buckets = {}
lines = 0
badlines = 0
with open(CLOSER_FILE) as closer_file:
	with open(GENE_LENGTHS_FILE,"w") as outfile:
		for line in closer_file:
			prb,cpg,distance,status = line.upper().split(",")[:4]
			if not is_missing(cpg,cpg_positions,missingcpg,"Missing cpg"):
				cpgchrm,cpgpos = cpg_positions[cpg]
				if not is_missing(prb,probemap,missingprobe,"Missing corresponding gene for probe"):
					gene = probemap[prb]
					if not is_missing(gene,genenames,missinggene,"Missing gene info for gene"):
						if (cpg,gene) not in sig_pairs:
							sig_pairs[(cpg,gene)] = status
			lines += 1

#now we need to load all this into a database for further analysis
conn = sqlite3.connect(":memory:")
c = conn.cursor()
c.execute("CREATE TABLE genes (cpg text,gene text)")
c.executemany("INSERT INTO genes VALUES (?,?)",sig_pairs)
c.execute("CREATE INDEX gene_name ON genes (gene)")
c.execute("CREATE INDEX cpg_name ON genes (cpg)")
conn.commit()

def query(stmt,fields):
	c = conn.cursor()
	c.execute(stmt,[fields])
	results = [x[0] for x in c.fetchall()]
	return results

def get_assoc_cpgs(gene):
	return query("SELECT DISTINCT cpg FROM genes WHERE gene = ?",gene)

def get_assoc_genes(cpg):
	return query("SELECT DISTINCT gene FROM genes WHERE cpg = ?",cpg)

c = conn.cursor()
stmt = "SELECT DISTINCT gene FROM genes" 
c.execute(stmt)	
genes = [x[0] for x in c.fetchall()]

groupid = 0
nodeid = 0

def getcpgbyname(name,cpgs):
	for cpg in cpgs:
		if "name" in cpg and cpg["name"] == name:
			return cpg
	return None

def getcpgtype(cpg,gene):
	global sig_pairs	
	status = sig_pairs[(cpg,gene)]
	cpgchrm,cpgpos = cpg_positions[cpg]
	genestart,genestop,genestrand,genechrm = genenames[gene]
	if cpgchrm != genechrm[3:]:
		return "trans"
	else:
		if (genestart >= cpgpos and genestop <= cpgpos) or (genestart <= cpgpos and genestop >= cpgpos):
			dist = 0
		else:
			dist = min(abs(genestart - cpgpos),abs(genestop-cpgpos))

		if dist >= 50000:
			return "distal"

		if dist == 0 or sig_pairs[(cpg,gene)] == "CLOSEST":
			return "in"

		return "notclosest"


def size(items):
	if not items or len(items) <= 0:
		return 1
	elif len(items) < 20:
		return len(items)
	elif len(items) < 100:
		return int(len(items) / 20) + 20
	else:
		return int(log(len(items))) + 30


#node_fname = open("GTP_Nodes.txt","w")
#link_fname = open("GTP_Links.txt","w")
node_entries = []
link_entries = []
def addgene(gene,groupid,nodes,links,cpgs=set()):
	global nodeid
	nodecnt = 0
	#dont process a cpg more than one time
	if gene in genes:
		genes.remove(gene)
		newcpgs = get_assoc_cpgs(gene)
		#create a node for the cpg, as we know its not in the list of nodes yet (otherwise we could ahve removed it from the list of cpgs)
		geneid = nodeid
		#node_fname.write(",".join([cpg,str(groupid),str(nodeid),str(len(genes))]) + "\n")
		#add the gene to the network
		node = {"id":gene,"node":nodeid,"type":"gene","group":groupid,"value":size(newcpgs)}
		nodes.append(node)
		nodeid += 1

		#now go through all the genes
		for cpg in newcpgs:
			cpggenes = get_assoc_genes(cpg)
			#get the node for the gene if its in the table, otherwise create it and increment the nodeid
			if cpg not in cpgs:
				nodecnt += 1
				cpgs.add(cpg)
			node = getcpgbyname(cpg,nodes)
			if not node:
				node  = getcpgbyname(cpg,node_entries)
				if not node:
					if len(cpggenes) == 1:
						cpgtype = getcpgtype(cpg,gene)
					else:
						cpgtype = "multi"
					node = {"id":"","type":cpgtype,"name":cpg,"node":nodeid,"group":groupid,"value":size(cpggenes)}
					nodes.append(node)
					nodeid += 1
			cpgid = node["node"]
			#now add a link between the cpg and the gene
			#link_fname.write(",".join([str(cpgid),str(geneid),str((len(cpgs) + len(genes)))]) + "\n")
			links.append({"source":cpgid,"target":geneid,"value":(size(newcpgs) + size(cpggenes)),"group":groupid,"cpg":cpg,"gene":gene})
			for nextgene in cpggenes:
				cnt = addgene(nextgene,groupid,nodes,links,cpgs)
				nodecnt += cnt
	return nodecnt

while len(genes) > 0:
	startid = nodeid
	nodes = []
	links = []
	nodecnt = addgene(genes[0],groupid,nodes,links,set())
	if nodecnt > 10:
		if groupid != 4:
			print "adding group %d of size %d" % (groupid, nodecnt)
			node_entries = node_entries + nodes
			link_entries = link_entries + links
		else:
			nodeid = startid
		groupid += 1
	else:
		nodeid = startid

with open("GTP_graph.json","w") as outfile:
	json.dump({"nodes":node_entries,"links":link_entries},outfile)

#now we need to use the links (supplemented with other information) to create circos png for all the groups that we've saved
for group in range(0,groupid):
	#get all the links that we have created for this particular group
	links = [x for x in link_entries if x["group"] == group]
	#now for each link, add the supplemental info to the data file that circos will read
	with open("links.txt","w") as outfile:
		with open("genes.txt","w") as genefile:
			genes = set()
			#chrm,start,stop,chrm,start,stop (of target and destination, space seperated)
			for link in links:
				cpgchrm,cpgpos = cpg_positions[link["cpg"]]
				genestart,genestop,genestrand,genechrm = genenames[link["gene"]]
				outfile.write(" ".join([str(x) for x in ["hs"+cpgchrm, cpgpos,int(cpgpos)+1,formatchrome(genechrm),genestart,genestop]]) + "\n")
				if link["gene"] not in genes:
						genes.add(link["gene"])
						genefile.write("%s\n" % (" ".join([str(x) for x in [formatchrome(genechrm),genestart,genestop,link["gene"]]])))
	
	#now that we've provided circos all teh info it needs, create the plot
	subprocess.call(["perl","/Users/georgegoehring/circos/circos-0.69/bin/circos","-conf","circos.conf"])
	#temporarily save the links for troubleshooting
	subprocess.call(["mv","links.txt","links" + str(group) + ".txt"])
	#rename the plot output to the correct filename for the group
	subprocess.call(["mv","circos.png",str(group) + ".png"])
	
