#!/usr/bin/python3

import sys, subprocess, os, re, csv, random, copy, argparse, atexit
from operator import itemgetter

def main(argv):

    args = get_args(argv)

    make_dirs()

    for entry in os.listdir(args.domain_dir):
        if not os.path.isdir(args.domain_dir + os.path.sep + entry):
            run_seqgen(args.domain_dir, entry, "domains", args.domain_length, args)

    for entry in os.listdir(args.gene_dir):
        if not os.path.isdir(args.gene_dir + os.path.sep + entry):
            run_seqgen(args.gene_dir, entry, "genes", args.gene_length, args)

    append_seqs("seqs" + os.path.sep + "domains", "seqs" + os.path.sep + "genes", args.leafmap_dir, args.append_domain)

#-----------------------------------------------------------------------------#

def get_args(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("domain_dir", help="Directory of input domain trees.")
    parser.add_argument("gene_dir", help="Directory of input gene trees.")
    parser.add_argument("leafmap_dir", help="Directory of domain to gene leaf maps.")
    parser.add_argument("-dl", "--domain-length", default="100", help="Length of domain sequences [default = 100].")
    parser.add_argument("-gl", "--gene-length", default="1000", help="Length of gene sequences [default = 1000].")
    parser.add_argument("-ap", "--append-domain", help="Append domain sequences to gene sequences instead of random insertion.", action="store_true")
    parser.add_argument("-s", "--branch-scaling", default="1.0", help="Branch length scaling factor [default = 1.0].")
    parser.add_argument("-m", "--model", default="GTR", help="HKY, F84, GTR, JTT, WAG, PAM, BLOSUM, MTREV, CPREV45, MTART, LG, and GENERAL. HKY, F84 and GTR are for nucleotides, while the rest are for amino acids [default = GTR].")
    parser.add_argument("-a", "--alpha", default="1.0", help="Shape (alpha) for gamma rate heterogeneity [default = 1.0].")
    parser.add_argument("-g", "--gamma-cats", help="Number of gamma rate categories [default = continuous].")
    parser.add_argument("-i", "--invariable-sites", default="0.0", help="Proportion of invariable sites [default = 0.0].")
    parser.add_argument("-c", "--codon", help="#1 #2 #3 = Rates for codon position heterogeneity [default = none].")
    parser.add_argument("-t", "--transition-transversion", default="1.0", help="Transition-transversion ratio [default = equal rate].")
    parser.add_argument("-r", "--rate-matrix", help="#1 #2 #3 #4 #5 #6 = General rate matrix [default = all 1.0].")
    parser.add_argument("-f", "--char-frequencies", default="e", help="#A #C #G #T = Nucleotide frequencies [default = all equal] or #1 .. #20 = Amino acid frequencies e = equal [default = matrix freqs].")
    parser.add_argument("-z", "--seed", help="Seed for random number generator [default = system generated].")
    args = parser.parse_args()
    return args

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

def run_seqgen(treedir, tree, dirtype, seqlen, args):

    infile_old = open(treedir + os.path.sep + tree, "r")
    line = infile_old.readlines()[0]
    line = re.sub(r'\[.*?\]', '', line.strip())
    line = re.sub(r'(\)[^:]*:)', r'):', line.strip())

    with open(treedir + os.path.sep + tree + "_temp", "w") as infile_new:
        infile_new.write(line.strip())

    infile = open(treedir + os.path.sep + tree + "_temp", "r")

    outfile = open("seqs" + os.path.sep + dirtype + os.path.sep + (("pre_transfer" + os.path.sep) if dirtype == "genes" else ("")) + tree, "w")

    call_list = ["./seq-gen", "-l", seqlen, "-s", args.branch_scaling, "-m", args.model, "-a", args.alpha, "-i", args.invariable_sites, "-of"]

    if args.gamma_cats != None:
        call_list.append("-g")
        call_list.append(args.gamma_cats)

    if args.codon != None:
        call_list.append("-c")
        call_list.append(args.codon)

    if args.model == "HKY" or args.model == "F84":
        call_list.append("-t")
        call_list.append(args.transition_transversion)

    if args.rate_matrix != None:
        call_list.append("-r")
        call_list.append(args.rate_matrix)

    if args.seed != None:
        call_list.append("-z")
        call_list.append(args.seed)

    subprocess.call(call_list, stdin=infile, stdout=outfile)

    os.remove(treedir + os.path.sep + tree + "_temp")

#-----------------------------------------------------------------------------#

def append_seqs(domainseqs, geneseqs, leafmaps, append_domain):

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
                #print("Trying insertion in Tree: " + str(col[2]) + " and Gene: " + str(col[1]) + " with History: " + str(genedict[col[2]][col[1]][1]))
                domain = domaindict[entry.split("leafmap")[0].strip() + "tree"][col[0]]
                gene = genedict[col[2]][col[1]][0]
                if append_domain:
                    genedict[col[2]][col[1]] = (genedict[col[2]][col[1]][0] + domain, genedict[col[2]][col[1]][1])
                else:
                    gene_adds = genedict[col[2]][col[1]][1]
                    insert_at = (random.randint(0, len(gene) + 1), len(domain))
                    sorted_adds = copy.deepcopy(gene_adds)
                    sorted_adds.append(insert_at)
                    sorted_adds.sort(key=itemgetter(0))
                    insert_index = sorted_adds.index(insert_at)
                    while (insert_index != 0 and (sorted_adds[insert_index - 1][0] + sorted_adds[insert_index - 1][1] - 1) >= insert_at[0]):
                        if insert_at[0] == 0 or insert_at[0] == len(gene):
                            break
                        #print(insert_at[0])
                        insert_at = (random.randint(0, len(gene) + 1), len(domain))
                        sorted_adds = copy.deepcopy(gene_adds)
                        sorted_adds.append(insert_at)
                        sorted_adds.sort(key=itemgetter(0))
                        insert_index = sorted_adds.index(insert_at)
                    #print("Insertion at: " + str(insert_at[0]) + " in Tree: " + str(col[2]) + " and Gene: " + str(col[1]))
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
