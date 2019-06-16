"""CopyNum - outputs the copy numbers for a list of input genes.

Usage:
  copyNum.py --input=<gene_list>
  copyNum.py --input <gene_list> --output <copy_numbers> --log <log_file>
  copyNum.py --input <gene_list> --seg <segment_file> --output <copy_numbers> --log <log_file>
  copyNum.py --input <gene_list> --genedef <gene_def> --seg <segment_file> --output <copy_numbers> --log <log_file>

Options:
  -h --help                Show this screen.
  --input=<in>             Filename of file containing a column of gene names.
  --output=<cn>            Filename of output file [default: copyNum.out.txt].
  --genedef=<gd>           Filename of gene file containing gene ranges [default: CCDS.current.txt].
  --seg=<sf>               Filename of seg file containing copy number data [default: CCLE_copynumber_2013-12-03.seg.txt].
  --log=<log>              Filename of log file [default: copyNum.log].
"""

"""
 * Author: Joshua Loving <joshua.k.loving@gmail.com>
 * Date:   February 16, 2019
 * Version: 0.1 (alpha)
"""

from docopt import docopt
import csv
from intervaltree import IntervalTree, Interval

def readInput(infile,genefile, segfile):
    """
    Reads input files.

    Extended description of function.

    Parameters:
    infile (str):   File containing list of genes to be analyzed
    genefile (str): File containing gene range definitions
    segfile (str):  File containing cell line intervals and copy number data

    Returns:
    genes (list):        List of genes
    genedef (dict):      Dictionary of genes mapping to corresponding intervals
    interval_dict(dict): Dictionary of dictionary of interval trees containing cell line ranges
    """
    with open(infile) as inf:
        genes = [i.strip() for i in inf.readlines()]
    with open(genefile) as genef:
        dictgenes = csv.DictReader(genef, delimiter="\t")
        genedef = {}
        for d in dictgenes:
            if d["cds_from"] != "-" and d["cds_to"] != "-":
                genedef[d["gene"]]  = (d["#chromosome"], Interval(int(d["cds_from"]),int(d["cds_to"])))

    with open(segfile) as seg:
        interval_dict = {}
        dictseg = csv.DictReader(seg, delimiter="\t")
        for d in dictseg:
            d = dict(d)
            if "e" in d["End"]:
                #Replace one incorrect exponential value
                d["End"] = 115000000
            if d["CCLE_name"] in interval_dict:
                if d["Chromosome"] in interval_dict[d["CCLE_name"]]:
                    interval_dict[d["CCLE_name"]][d["Chromosome"]][int(d["Start"]):int(d["End"])] = float(d["Segment_Mean"])
                else:
                    interval_dict[d["CCLE_name"]][d["Chromosome"]] = IntervalTree()
                    interval_dict[d["CCLE_name"]][d["Chromosome"]][int(d["Start"]):int(d["End"])] = float(d["Segment_Mean"])
            else:
                interval_dict[d["CCLE_name"]] = dict()
                interval_dict[d["CCLE_name"]][d["Chromosome"]] = IntervalTree()
                interval_dict[d["CCLE_name"]][d["Chromosome"]][int(d["Start"]):int(d["End"])] = float(d["Segment_Mean"])

    return genes, genedef, interval_dict

import datetime

if __name__ == '__main__':
    #Parse arguments from the docstring
    arguments = docopt(__doc__)

    #Generate logfile
    with(open(arguments["--log"], "w")) as logf:
        logf.write("Date: %s\n"% str(datetime.date.today()))
        logf.write("Version: CopyNum 0.1")
        logf.write("Input gene file: %s\n" % arguments["--input"])
        logf.write("Output CN file: %s\n" % arguments["--output"])
        logf.write("Input gene range file: %s\n" % arguments["--genedef"])
        logf.write("Input cell line CN file: %s\n" % arguments["--seg"])

    #parse input into objects
    g, gdef, id = readInput(arguments["--input"],arguments["--genedef"],arguments["--seg"])

    #Get only the genes we are interested in
    interest = dict((k, gdef[k]) for k in g if k in gdef)

    with open(arguments["--output"], "w") as out:
        #Generate header
        out.write("Gene\t")
        for line in id:
            out.write(line + "\t")
        out.write("\n")
        #Loop over genes
        for gene_name in interest:
            gene = interest[gene_name]

            out.write(gene_name + "\t")
            #Loop over cell line intervals
            for line in id:

                SumCN = 0
                SumWidth = 0

                #Check whether there were overlapping intervals found
                if len(id[line][gene[0]][gene[1].begin:gene[1].end]) == 0:
                    #No copy number data, write NA
                    out.write("NA\t")
                else:
                    for i in id[line][gene[0]][gene[1].begin:gene[1].end]:
                        #Compute copy number across interval
                        width = i.end - i.begin
                        SumCN += i.data*width
                        SumWidth += width
                    #Average copy number across gene intervals
                    out.write(str(SumCN/SumWidth))
                    out.write("\t")
            out.write("\n")
