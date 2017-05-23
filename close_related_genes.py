import sqlite3

#this should take in a closer file, load it into sql, and then query it for genes that are related to the same cpg and overlapping or close to each other

CLOSER_GENES_FILE = "/home/ekennedy/Projects/Thesis/output/GTP_closer" + ".csv"
MAX_DISTANCE = 2000

def load_relations(fname):
	with open(fname) as infile:
		for line in infile:
			exp_probe,cpg,chrm,distance,status,gene,strand,gene_start,gene_end,gene_dist,pval,tstat,beta,beta_sd = line.strip().split(",")
			if strand == '-':
				tmp = gene_end
				gene_end = gene_start
				gene_start = tmp
			yield (exp_probe,cpg,chrm,int(distance),status,gene,strand,int(gene_start),int(gene_end),int(gene_dist),float(pval),float(tstat),float(beta),float(beta_sd))
	
conn = sqlite3.connect(":memory:")
c = conn.cursor()
c.execute("CREATE TABLE genes (exp_probe text,cpg text,chrm text, distance int, status text,gene text,strand text,start int, end int, gene_dist int,pval float,tstat float,beta float,beta_sd float)")
c.executemany("INSERT INTO genes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",load_relations(CLOSER_GENES_FILE))
c.execute("CREATE INDEX	chrm_cpg ON genes (chrm,cpg)")
c.commit()
		
overlap = [
"g.strand = '+' AND o.strand = '+' AND g.start <= o.start AND g.end > o.start"
#
#  ==G.START==O.START==O.END==G.END==>
#  <==================================
#
# AND
#
#  ==G.START==O.START==G.END==O.END==>
#  <==================================
#
#  OVERLAP: MIN(O.END,G.END) - O.START

,"g.strand = '-' AND o.strand = '-' AND g.start >= o.start AND g.end <= o.start"
#
#  ==================================>
#  <=G.END==O.END==O.START==G.START===
#
#  AND
#
#  ==================================>
#  <=O.END==G.END==O.START==G.START==
#
#  OVERLAP: O.START - MAX(G.END,O.END)


,"g.strand = '-' AND o.strand = '+' AND g.start >= o.end AND g.end <= o.end"
#
#  =========O.START==O.END==========>
#  <=G.END==================G.START==
#
#  AND
#
#  ==O.START=========O.END===========>
#  <==========G.END=========.G.START==
#
#  OVERLAP: O.END - MAX(O.START,G.END)

,"g.strand = '+' AND o.strand = '-' AND g.start <= o.start AND g.end >= o.start"
#
#  ==G.START==================G.END==>
#  <==========O.END==O.START==========
#
#  AND
#
#  =========G.START===========G.END==>
#  <=O.END===========O.START==========
#
#  OVERLAP: O.START - MAX(G.START,O.END)
]

distances = [
"g.strand = '+' AND o.strand = '+' AND g.end < o.start AND o.start-g.end < %d" %(MAX_DISTANCE)
#
#  ==G.START==G.END==O.START==O.END==>
#  <==================================

,"g.strand = '-' AND o.strand = '-' AND g.end > o.start AND g.end - o.start < %d" %(MAX_DISTANCE)
#
#  ==================================>
#  <=O.END==O.START==G.END==G.START===

,"g.strand = '+' AND o.strand = '-' AND g.end < o.end AND o.end - g.end < %d" %(MAX_DISTANCE)
#
#  ==G.START==G.END==================>
#  <=================O.END==O.START===

,"g.strand = '+' AND o.strand = '-' AND o.start < g.start AND g.start - o.strand < %d" %(MAX_DISTANCE)
#
#  ==================G.START==G.END==>
#  <=O.END==O.START===================
]

def cmbn_stmts(stmts):
	return " OR ".join(["("+x+")" for x in stmts])

query = "SELECT g.exp_probe,g.cpg,g.chrm,g.distance,g.status,g.gene,g.strand,g.gene_start,g.gene_end,g.gene_dist,g.pval,g.tstat,g.beta,g.beta_sd FROM genes g,genes o WHERE g.chrm = o.chrm AND g.cpg = o.cpg AND (%s)" % cmbn_stmts(overlap) # + distances)

c = conn.cursor()
c.executequery(query)
with open(OUTFILE,'w') as outfile:
	for result in c.fetchall():
		line = "%s\n" % (",".join([str(x) for x in result]))
		outfile.write(line)
