import sys
cols = {}
for line in sys.stdin:
	cnt = len(line.split("\t"))
	cols[cnt] = cols.get(cnt,0) + 1

for k,v in cols.iteritems():
	print "%d\t%d" % (k,v)
