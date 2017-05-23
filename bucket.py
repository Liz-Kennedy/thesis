import sys
buckets = {}
for line in sys.stdin:
	buckets[line] = buckets.get(line,0) + 1
for k,v in buckets.iteritems():
	print k.strip() + "," + str(v)
