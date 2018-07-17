#!/usr/bin/python3

import sys, subprocess, os, re, csv, random
from operator import itemgetter

def main(argv):

    parse_args(argv)

    make_dirs()

    for entry in os.listdir(argv[0]):
        if not os.path.isdir(argv[0] + os.path.sep + entry):
            run_seqgen(argv[0], entry, "domains", "100")

    for entry in os.listdir(argv[1]):
        if not os.path.isdir(argv[1] + os.path.sep + entry):
            run_seqgen(argv[1], entry, "genes", "1000")

    append_seqs("seqs" + os.path.sep + "domains", "seqs" + os.path.sep + "genes", argv[2])

#-----------------------------------------------------------------------------#

def parse_args(argv):
    pass

#-----------------------------------------------------------------------------#

def make_dirs():

    if not os.path.exists("seqs"):
        os.mkdir("seqs")

    if not os.path.exists("seqs" + os.path.sep + "domains"):
        os.mkdir("seqs" + os.path.sep + "domains")

    if not os.path.exists("seqs" + os.path.sep + "genes"):
        os.mkdir("seqs" + os.path.sep + "genes")

    if not os.path.exists("seqs" + os.path.sep + "genes" + os.path.sep + "pre_transfer"):
        os.mkdir("seqs" + os.path.sep + "genes" + os.path.sep + "pre_transfer")

    if not os.path.exists("seqs" + os.path.sep + "genes" + os.path.sep + "post_transfer"):
        os.mkdir("seqs" + os.path.sep + "genes" + os.path.sep + "post_transfer")

#-----------------------------------------------------------------------------#

def run_seqgen(treedir, tree, dirtype, seqlen):

    infile_old = open(treedir + os.path.sep + tree, "r")
    line = infile_old.readlines()[0]
    line = re.sub(r'\[.*?\]', '', line.strip())
    line = re.sub(r'(\)[^:]*:)', r'):', line.strip())

    with open(treedir + os.path.sep + tree + "_temp", "w") as infile_new:
        infile_new.write(line.strip())

    infile = open(treedir + os.path.sep + tree + "_temp", "r")

    outfile = open("seqs" + os.path.sep + dirtype + os.path.sep + (("pre_transfer" + os.path.sep) if dirtype == "genes" else ("")) + tree, "w")
    subprocess.call(["../seq-gen", "-l", seqlen, "-m", "GTR", "-a", "1", "-of"],
                    stdin=infile, stdout=outfile)

    os.remove(treedir + os.path.sep + tree + "_temp")

#-----------------------------------------------------------------------------#

def append_seqs(domainseqs, geneseqs, leafmaps):

    domaindict = {}
    genedict = {}

    for entry in os.listdir(domainseqs):
        domaindict[entry] = {}
        latest_key = ""
        with open(domainseqs + os.path.sep + entry, "r") as domainfile:
            for line in domainfile:
                if line.startswith(">"):
                    latest_key = line.split(">")[1].strip()
                    domaindict[entry][latest_key] = ""
                else:
                    domaindict[entry][latest_key] += line.strip()

    for entry in os.listdir(geneseqs + os.path.sep + "pre_transfer"):
        genedict[entry] = {}
        latest_key = ""
        with open(geneseqs + os.path.sep + "pre_transfer" + os.path.sep + entry, "r") as genefile:
            for line in genefile:
                if line.startswith(">"):
                    latest_key = line.split(">")[1].strip()
                    genedict[entry][latest_key] = ("",[])
                else:
                    genedict[entry][latest_key] = (genedict[entry][latest_key][0] + line.strip(), [])

    for entry in os.listdir(leafmaps):
        with open(leafmaps + os.path.sep + entry, "r") as leafmapfile:
            for col in csv.reader(leafmapfile, delimiter="\t"):
                domain = domaindict[entry.split("leafmap")[0].strip() + "tree"][col[0]]
                gene = genedict[col[2]][col[1]][0]
                gene_adds = genedict[col[2]][col[1]][1]
                insert_at = (random.randint(0, len(gene) + 1), len(domain))
                sorted_adds = gene_adds
                sorted_adds.append(insert_at)
                sorted_adds.sort(key=itemgetter(0))
                insert_index = sorted_adds.index(insert_at)
                while ((sorted_adds[insert_index - 1][0] + sorted_adds[insert_index - 1][1] - 1) >= insert_index):
                    if insert_at[0] == 0 or insert_at[0] == len(gene):
                        break
                    insert_at = (random.randint(0, len(gene) + 2 * len(domain)), len(domain))
                    sorted_adds = gene_adds
                    sorted_adds.append(insert_at)
                    sorted_adds.sort(key=itemgetter(0))
                    insert_index = sorted_adds.index(insert_at)
                genedict[col[2]][col[1]] = (genedict[col[2]][col[1]][0][:insert_at[0]] + domain + genedict[col[2]][col[1]][0][insert_at[0]:], sorted_adds)
                new_tuples = genedict[col[2]][col[1]][1]
                for i,j in enumerate(new_tuples):
                    if i > insert_index:
                        new_tuples[i] = (j[0] + 1, j[1])
                genedict[col[2]][col[1]] = (genedict[col[2]][col[1]][0], new_tuples)

    for entry in os.listdir(geneseqs + os.path.sep + "pre_transfer"):
        with open(geneseqs + os.path.sep + "post_transfer" + os.path.sep + entry, "w") as newgenefile:
            for key in genedict[entry]:
                newgenefile.write(">" + key + "\n")
                newgenefile.write(genedict[entry][key][0] + "\n")

#-----------------------------------------------------------------------------#

if __name__ == "__main__":
    main(sys.argv[1:])
