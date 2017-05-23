import sqlite3
from random import randint
DEBUG = False
DATASET="GTP"# or "MESA"

def loadprobelocs(filename):
        probes = []
        with open(filename) as infile:
                for line in infile:
                        #IlluminaID, CHRM, Probe.start, Probe.end, START, END, STRAND, pref_target, GENE
			nm,chrm,prb_start,prb_end,start,stop,strand,category,gene = line.strip().split("\t")
                        #if the strand is negative we need to swap the start and stop locations
			if strand == "-":
                                tmp = stop
                                stop = start
                                start = tmp
                        if start != "NA" and stop != "NA":
				if gene.upper().startswith("ILMN"):
					annotated = 1
				else:
					annotated = 0
				probes.append((nm,frmt_chrm(chrm),int(prb_start),int(prb_end),int(start),int(stop),strand,format_category(category),gene,annotated))
        load(probes)

def format_category(cat):
	cat = cat.upper()
	if "PREF" in cat:
		return 0
	elif "SEC" in cat:
		return 1
	else:
		return 2

def closest_probe(probe,chrm,pos):
	chrm = frmt_chrm(chrm)
	sql = """
	SELECT chrm,prb_start,prb_end,strand,gene,start,end,start_dist,category,
		CASE WHEN chrm IS null THEN 4
		WHEN start_dist = 0 THEN 0
		WHEN start_dist < 2500 THEN 1
		WHEN start_dist < 1000000 THEN 2
		ELSE 3 END as priority
	FROM 
	(SELECT chrm,probe,gene,prb_start,prb_end,strand,start,end,category,annotated,
		CASE WHEN chrm <> ? THEN null
	 	WHEN (strand = '+' AND start <= ? AND end >= ?) OR (strand = '-' AND start >= ? AND end <= ?) THEN 0
		ELSE ABS(start - ?) END as start_dist
		FROM probes WHERE probe = ?) tmp
	ORDER BY priority,annotated,category,start_dist ASC
	"""
	c = conn.cursor()
	c.execute(sql,(chrm,pos,pos,pos,pos,pos,probe))
	results = c.fetchall()
	if DEBUG:
		print "---RESULTS---"
		for result in results:
			result = [str(x) for x in result]
			print ",".join(result)
		print "-------------"
	if len(results) == 0:
		return None #raise Exception("No matching probes for %s" % (probe))#LIZ: I have removed some probes from analysis. they have no match.**********
	else:
		return results[0]

def frmt_chrm(chrm):
	chrm = chrm.lower()
	if chrm.startswith("chr"):
		return chrm[3:]
	else:
		return chrm	
conn = None

def load(probes):
	global conn 
	conn = sqlite3.connect(":memory:")
	c = conn.cursor()
	c.execute("CREATE TABLE probes (probe text,chrm text, prb_start int, prb_end int, start int, end, strand varchar(1), category int,gene text, annotated int)")
	c.executemany("INSERT INTO probes VALUES (?,?,?,?,?,?,?,?,?,?)",probes)
	c.execute("CREATE INDEX probe_name ON probes (probe)")
	c.execute("CREATE INDEX probe_loc ON probes (chrm,start,end)")
	c.execute("CREATE INDEX probes_strand ON probes (chrm,strand)")
	conn.commit()

def drop():
	if conn:
		c = conn.cursor()
		c.execute("DROP TABLE IF EXISTS probes")
		conn.commit()
