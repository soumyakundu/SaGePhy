#!/usr/bin/python3

import sys, subprocess, os, re, csv

def main(argv):

    parse_args(argv)
    run_seqgen(os.getcwd(), argv[0])

    if not os.path.exists(argv[1] + os.path.sep + "seqs"):
        os.mkdir(argv[1] + os.path.sep + "seqs")

    for entry in os.listdir(argv[1]):
        if not os.path.isdir(argv[1] + os.path.sep + entry):
            run_seqgen(argv[1], entry)

    append_seqs(argv[0] + "_seqs", argv[1] + os.path.sep + "seqs", argv[2])

#-----------------------------------------------------------------------------#

def parse_args(argv):
    pass

#-----------------------------------------------------------------------------#

def run_seqgen(treedir, tree):

    infile_old = open(treedir + os.path.sep + tree, "r")
    line = infile_old.readlines()[0]
    line = re.sub(r'\[.*?\]', '', line.strip())
    line = re.sub(r'(\)[^:]*:)', r'):', line.strip())

    with open(treedir + os.path.sep + tree + "_temp", "w") as infile_new:
        infile_new.write(line.strip())

    infile = open(treedir + os.path.sep + tree + "_temp", "r")

    if treedir != os.getcwd():
        outfile = open(treedir + os.path.sep + "seqs" + os.path.sep + tree + "_seqs", "w")
        subprocess.call(["../seq-gen", "-m", "GTR", "-a", "1", "-of"],
                     stdin=infile, stdout=outfile)
    else:
        outfile = open(tree + "_seqs", "w")
        subprocess.call(["../seq-gen", "-l", "100", "-m", "GTR", "-a", "1", "-of"],
                        stdin=infile, stdout=outfile)

    os.remove(treedir + os.path.sep + tree + "_temp")

#-----------------------------------------------------------------------------#

def append_seqs(domainseqs, geneseqsdir, leafmap):

    domaindict = {}
    genedict = {}

    latest_key = ""
    with open(domainseqs, "r") as domainfile:
        for line in domainfile:
            if line.startswith(">"):
                latest_key = line.split(">")[1].strip()
                domaindict[latest_key] = ""
            else:
                domaindict[latest_key] += line.strip()

    for entry in os.listdir(geneseqsdir):
        genedict[entry] = {}
        latest_key = ""
        with open(geneseqsdir + os.path.sep + entry, "r") as genefile:
            for line in genefile:
                if line.startswith(">"):
                    latest_key = line.split(">")[1].strip()
                    genedict[entry][latest_key] = ""
                else:
                    genedict[entry][latest_key] += line.strip()

    with open(leafmap) as tsv:
        for col in csv.reader(tsv, delimiter="\t"):
            genedict[col[2] + "_seqs"][col[1]] += domaindict[col[0]]

    for entry in os.listdir(geneseqsdir):
        with open(geneseqsdir + os.path.sep + entry, "w") as newgenefile:
            for key in genedict[entry]:
                newgenefile.write(">" + key + "\n")
                newgenefile.write(genedict[entry][key] + "\n")

#-----------------------------------------------------------------------------#

if __name__ == "__main__":
    main(sys.argv[1:])
