import gffutils
import sys
from liftofftools.cli_arguments import ARGS


def get_database_connection(gff_file):
    try:
        gffutils.FeatureDB(gff_file)
    except:
        feature_db = build_database(gff_file)
    feature_db.execute('ANALYZE features')
    return feature_db


def build_database(gff_file):
    print(ARGS.infer_genes)
    try:
        feature_db = gffutils.create_db(gff_file, gff_file + "_db", merge_strategy="create_unique",
                                        force=True, verbose=True, checklines=10, disable_infer_genes=ARGS.infer_genes)
    except:
        find_problem_line(gff_file)
    return feature_db



def find_problem_line(gff_file):
    f = open(gff_file, 'r')
    lines = f.readlines()
    for i in range (len(lines)):
        line = lines[i]
        if line[0] != "#":
            try:
                gffutils.create_db(line, ":memory", from_string=True, force=True)
            except:
                sys.exit("ERROR:Incorrect GFF/GTF syntax on line " + str(i + 1))