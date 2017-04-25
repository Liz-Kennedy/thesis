import sqlite3

DATASET="GTP"# or "MESA"

def loadprobelocs(filename):
        probes = []
        with open(filename) as infile:
                for line in infile:
                        nm,chrm,start,stop,strand,category = line.strip().split("\t")
                        #if the strand is negative we need to swap the start and stop locations
                        if strand == "-":
                                tmp = stop
                                stop = start
                                start = tmp
                        if start != "NA" and stop != "NA":
				probes.append((nm,frmt_chrm(chrm),int(start),int(stop),strand,format_category(category)))
        load(probes)

def format_category(cat):
	cat = cat.upper()
	if "PRI" in cat:
		return 0
	elif "SEC" in cat:
		return 1
	else:
		return 2

def closest_probe(probe,chrm,pos):
	chrm = frmt_chrm(chrm)
	sql = """SELECT chrm,start,end,strand,(CASE WHEN chrm = ? THEN ABS(start - ?) ELSE null END) AS start_dist,category,
	(CASE WHEN chrm = ? THEN
		CASE WHEN start <= ? AND end >= ? THEN 0 ELSE
			CASE WHEN ((strand = '+' AND ? < start) OR (strand = '-' AND ? > start)) AND ABS(start - ?) < 1500 THEN 1 ELSE 2 END
		END
	ELSE 3 END) as priority
	FROM probes WHERE probe = ? ORDER BY priority,category,start_dist ASC
	"""
	c = conn.cursor()
	c.execute(sql,(chrm,pos,chrm,pos,pos,pos,pos,pos,probe))
	results = c.fetchall()
	if len(results) == 0:
		raise Exception("No matching probes for %s" % (probe))
	elif len(results) == 1 or (len(results) > 0 and results[0][6] < 3):
		return results[0]
	else:
		category = results[0][5]
		results = [x for x in results if x[5] == category]
		return results[randint(len(results))]

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
	c.execute("CREATE TABLE probes (probe text,chrm text, start int, end int, strand text, category int)")
	c.executemany("INSERT INTO probes VALUES (?,?,?,?,?,?)",probes)
	c.execute("CREATE INDEX probe_name ON probes (probe)")
	c.execute("CREATE INDEX probe_loc ON probes (chrm,start,end)")
	c.execute("CREATE INDEX probes_strand ON probes (chrm,strand)")
	conn.commit()

def drop():
	if conn:
		c = conn.cursor()
		c.execute("DROP TABLE IF EXISTS probes")
		conn.commit()
