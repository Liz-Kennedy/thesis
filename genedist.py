import sqlite3

DATASET="MESA"# or "GTP"

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

def dist_type(dist):
	if(dist == 0): return "IN"
	elif (dist < 0): return "UPSTREAM"
	else: return "DOWNSTREAM"

def frmt_chrm(chrm):
	chrm = chrm.lower()
	if chrm.startswith("chr"):
		return chrm[3:]
	else:
		return chrm	

conn = None
def load_file(fname):
	genes = []
	with open(fname) as infile:
		infile.readline() #drop the header
		for line in infile:
			if DATASET == "MESA":
				gene,chrm,start,end,strand = line.strip().split("\t")[:5]
			elif DATASET == "GTP":
				gene,chrm,strand,start,end = line.strip().split("\t")[:5]
			if start == "NA" or end == "NA" or chrm == "NA":
                                continue
			start = int(start)
			end = int(end)
			chrm = frmt_chrm(chrm)
			if strand == "-":
				tmp = start
				start = end
				end = tmp
			genes.append((gene,chrm,start,end,strand))
	load(genes)

def load(genes):
	global conn 
	conn = sqlite3.connect(":memory:")
	c = conn.cursor()
	c.execute("CREATE TABLE genes (gene text,chrm text, start int, end int, strand text)")
	c.executemany("INSERT INTO genes VALUES (?,?,?,?,?)",genes)
	c.execute("CREATE INDEX gene_name ON genes (gene)")
	c.execute("CREATE INDEX gene_loc ON genes (chrm,start,end)")
	c.execute("CREATE INDEX gene_strand ON genes (chrm,strand)")
	conn.commit()

def loadcpgs(cpgs):
	global conn
	if not conn:
		conn = sqlite3.connect(":memory:")
	c = conn.cursor()
	c.execute("CREATE TABLE cpgs (cpg text,chrm text,cpgpos int)")
	c.executemany("INSERT INTO cpgs VALUES (?,?,?)",cpgs)
	c.execute("CREATE INDEX cpg_posinfo ON cpgs (chrm,cpgpos)")
	conn.commit()

def all_genes():
	c = conn.cursor()
	stmt = "SELECT gene,chrm,start,end,strand FROM genes"
	c = conn.cursor()
	c.execute(stmt)
	return c.fetchall()

def drop():
	if conn:
		c = conn.cursor()
		c.execute("DROP TABLE IF EXISTS genes")
		conn.commit()

#the start chrm and start position of the gene, and the position of the cpg
def findtss(chrm,start,pos):
	tss_dist = abs(start - pos)
	#find starting locations that are less abs distance away from the specified position
	return query("chrm = ? AND ABS(start - ?) < ?",(chrm,pos,tss_dist))

#chrm and distance to cpg for significant gene, and the position fo the cpg
def findcloser(chrm,dist,pos,limit=False):
	stmt = "chrm = ? AND ((start <= ? AND end >= ?) OR (end <= ? AND start >= ?) OR MIN(ABS(start - ?),ABS(end - ?)) < ?)"
	if limit:
		stmt += " LIMIT 1"
	return query(stmt,(chrm,pos,pos,pos,pos,pos,pos,dist))

#pass in a genes chrm,start, and distance to find genes that are between it and the cpg location
def findbetween(chrm,start,pos):
	if start < pos:
		return query("chrm = ? AND start > ? AND start < ?", (chrm,start,pos))
	else:
		return query("chrm = ? AND start > ? AND start < ?", (chrm,pos,start))

def find_abs_closest(cpgname,chrm,pos):
	c = conn.cursor()
	
	stmt = "SELECT dist FROM (SELECT CASE WHEN (start <= ? AND end >= ?) OR (start >= ? AND end <= ?) THEN 0 ELSE MIN(ABS(start-?),ABS(end-?)) END AS dist FROM genes WHERE chrm = ?) as tmp ORDER BY dist ASC LIMIT 1"

	return c.execute(stmt,[pos,pos,pos,pos,pos,pos,chrm]).fetchone()

def find_close_genes():
	c = conn.cursor()
	stmt = "SELECT gene,start,end,strand,cpg,cpgpos,chrm,dist FROM (SELECT gene,start,end,strand,cpg,cpgpos,cpgs.chrm as chrm,(CASE WHEN (start <= cpgpos AND end >= cpgpos) OR (start >= cpgpos AND end <= cpgpos) THEN 0 ELSE MIN(ABS(start-cpgpos),ABS(end-cpgpos)) END) as dist FROM genes,cpgs WHERE genes.chrm = cpgs.chrm) as distances WHERE ABS(dist) <= 2000 ORDER BY cpg,dist ASC"
	return c.execute(stmt).fetchall()

def within_range(chrm,pos,maxdist,mindist):
	c = conn.cursor()
	if mindist < 0:
		results = query("chrm = ? AND ((strand = '+' AND ?-start <= ? AND ?-start >= ?) OR (strand = '-' AND start - ? <= ? AND start - ? >= ?))",[chrm,pos,mindist,pos,maxdist,pos,mindist,pos,maxdist])
	else:
		results = query("chrm = ? AND ((strand = '+' AND ? - end >= ? AND ?-end <= ?) OR (strand = '-' AND end - ? >= ? AND end - ? <= ?))",[chrm,pos,mindist,pos,maxdist,pos,mindist,pos,maxdist]) 
	return len(results) > 0

def find_closest_tss(cpgname,chrm,pos):
	c = conn.cursor()
	
	stmt="SELECT dist FROM (SELECT ABS(start-?) AS dist FROM genes WHERE chrm = ?) as tmp ORDER BY dist ASC limit 1"

	return c.execute(stmt,[pos,chrm]).fetchone()

def query(where_clause,fields):
	c = conn.cursor()
	stmt = "SELECT gene,start,end,strand FROM genes WHERE %s" % (where_clause)
	c.execute(stmt,fields)
	return c.fetchall()

def stmt(query):
	c = conn.cursor()
	c.execute(query)
	return c.fetchall()
