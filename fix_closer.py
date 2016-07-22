import itertools
import sys

CPG_GENE_FILE = "MESA_closer.csv"
OUTFILE = "MESA_closer2.txt"
keyfn = lambda x: x[:3]

keep_status = set(["IN","CLOSEST","TRANS","DISTAL"]) 

PROBE_TO_GENE_MAPPING = "Illumina_HT12_genenames.txt"
probe_to_gene = {}
with open(PROBE_TO_GENE_MAPPING) as infile:
    for line in infile:
        prb,gene = line.strip().split(" ")
        probe_to_gene[prb] = gene

GENE_SELECTION_HISTORY = "fix_closer_selection_memory.txt"
prev_selections = {}
with open(GENE_SELECTION_HISTORY) as infile:
    for line in infile:
        gene,other_gene,selection = line.strip().split(",")
        prev_selections[(gene,other_gene)] = selection

missing = set()
def loadcpggenes():
        with open(CPG_GENE_FILE) as infile:
                for line in infile:
			#sig_gene,cpgname,distance,status,other_gene,other_distance,pval,fstat,beta,beta_sd
                        probe,cpg,distance,status,other_gene,other_distance,pval,fstat,beta,beta_sd = line.strip().split(",")
                        if probe in probe_to_gene:
				gene = probe_to_gene[probe]
                        	yield (cpg,probe,other_gene,other_distance,gene,distance,status,pval,fstat,beta,beta_sd)
			elif probe not in missing:
				missing.add(probe)
				print "WARN Missing corresponding gene for probe %s" % (probe)

#get the # of lines in the cpg gene file
cpg_gene_linecnt = len(list(loadcpggenes()))

def shouldkeep(entry):
    return "Y" #short circuits the brain
    cpg,probe,other_gene,other_distance,gene,distance,status = entry
    if (gene,other_gene) in prev_selections:
        line = prev_selections[(gene,other_gene)]
    else:
        print "\t".join([gene,distance,status,other_gene,other_distance])
	print "Keep [Y\N\X]:"
        line = sys.stdin.readline().strip().upper()
        if line == "X":
            exit()
        elif line != "Y":
            line = "N"
        prev_selections[(gene,other_gene)] = line
    return line == "Y"
    
def get_status(entries,status):
    entries = [x for x in entries if x[3] == status]
    if len(entries) == 0:
        return None
    else:
        return entries[0]

def getclosest(entries):
    for status in ["BETWEEN","TSS","CLOSER"]:
        match = get_status(entries,status)
        if match:
            if shouldkeep(match):
                return match
            else:
                return None
    return None

#save the programs brain before exiting
def exit():
    with open(GENE_SELECTION_HISTORY,"w") as outfile:
        for (gene,other_gene),response in prev_selections.iteritems():
            outfile.write("%s,%s,%s\n" % (gene,other_gene,response))
    sys.exit(0)

def getkeepers(entries):
    #filter out all the "other" genes that match the gene the probe is associated with
    diffgenes = [x for x in entries if x[2] != x[4]]
    #if the other distance is zero, drop the gene
    keepers = []
    for (cpg,prb,other_gene),entries in itertools.groupby(diffgenes,lambda x: x[:3]):
        entries = list(entries)
        closest = getclosest(entries)
        if closest:
            keepers.append(closest)
    return keepers

def writefn(entry,handle):
    cpg,probe,other_gene,other_distance,gene,distance,status,pval,fstat,beta,beta_sd = entry
    #gene,distance,other_gene,other_distance
    handle.write("\t".join([cpg,probe,gene,distance,status,other_gene,other_distance,pval,fstat,beta,beta_sd]) + "\n")

with open(OUTFILE,"w") as outfile:
        outfile.write("\t".join(["cpg","probe","gene","distance","status","other_gene","other_distance","pval","fstat","beta","beta_sd"]) + "\n")
	lines = 0
        #group the entries together that are associated with the same cpg & probe
        for (cpg,prb),entries in itertools.groupby(sorted(loadcpggenes()),lambda x: x[:2]):
           print "%s,%s %s percent complete" % (cpg,str(prb),str(100 * float(lines)/float(cpg_gene_linecnt)))
           entries = list(entries)
           lines += len(entries)
           #if the significant gene is the closest to the cpg, no checking, no 'other' genes are relevant, so write directly to output
           if entries[0][6] in keep_status:
               writefn(entries[0],outfile) 
           #if there were genes that fell closer
           else:
               keepers = getkeepers(entries)
               #if there were no keepers, then we can assume the significant gene is the closest
               if len(keepers) == 0: 
                   cpg,probe,other_gene,other_distance,gene,distance,status,pval,fstat,beta,beta_sd = entries[0]
                   writefn((cpg,probe,"NA","NA",gene,distance,"CLOSEST",pval,fstat,beta,beta_sd),outfile)
               else:
                   for keeper in keepers:
                       writefn(keeper,outfile)
exit()
