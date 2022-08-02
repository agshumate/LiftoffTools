from liftofftools import annotation
from liftofftools.sequences import SequenceDict
from liftofftools.cluster_analysis import analyze_clusters
from pyfaidx import Fasta
from liftofftools.synteny import synteny
from liftofftools.cli_arguments import parse_args
from liftofftools.filepaths import make_directory
from liftofftools.variants import variants
import sys
import warnings


def main(arglist=None):
    args= parse_args(arglist)
    make_directory(args.dir)
    feature_types = get_feature_types(args.f)
    args.ft = feature_types
    ref_db = annotation.Annotation(args.rg, args.infer_genes)
    target_db = annotation.Annotation(args.tg, args.infer_genes)
    ref_fa = Fasta(args.r)
    target_fa = Fasta(args.t)
    check_feature_types(target_db, feature_types)
    check_ids(ref_db.db_connection, target_db.db_connection, feature_types)
    child_types = get_child_types(feature_types, ref_db)
    if args.subcommand in ['all', 'synteny']:
        synteny.analyze_synteny(ref_db, target_db,ref_fa, target_fa, args)
    if args.subcommand in ['all', 'clusters', 'variants']:
        print('Extracting transcript sequences')
        ref_proteins = SequenceDict(ref_db, ref_fa, ['CDS', 'start_codon', 'stop_codon'], True)
        target_proteins = SequenceDict(target_db, target_fa, ['CDS','start_codon', 'stop_codon'], True)
        ref_trans = SequenceDict(ref_db, ref_fa, child_types, False)
        target_trans = SequenceDict(target_db, target_fa, child_types, False)
    if args.subcommand in ['all', 'clusters']:
        analyze_clusters.main(ref_proteins, target_proteins, ref_trans, target_trans, ref_db, target_db, args)
    if args.subcommand in ['all', 'variants']:
        variants.analyze_variants(ref_proteins, target_proteins, ref_trans, target_trans, target_db, ref_db, args)




def get_feature_types(feature_arg):
    feature_types = ['gene']
    if feature_arg is not None:
        with open(feature_arg) as fa:
            for line in fa:
                feature_types.append(line.strip())
    return feature_types


def get_child_types(parent_types, db):
    child_types = set()
    for parent in parent_types:
        for feature in db.db_connection.features_of_type(featuretype=parent):
            child_count = 0
            for child in db.db_connection.children(feature):
                child_count += 1
                if db.is_lowest_child(child):
                    child_types.add(child.featuretype)
            if child_count == 0:
                child_types.add(feature.featuretype)
    return child_types


def check_feature_types(db, feature_types):
    db_feature_types = list(db.db_connection.featuretypes())
    for f_type in feature_types:
        if f_type in db_feature_types:
            return
    sys.exit('No ' + ", ".join(feature_types) + " features found in reference annotation. Use -f option to "
                                                "provide text file of features types to analyze")


def check_ids(ref_db, target_db, feature_types):
    ref_ids = [feature.id for feature in ref_db.all_features() if feature.featuretype in feature_types]
    target_ids = [feature.id for feature in target_db.all_features() if feature.featuretype in feature_types]
    if set(ref_ids).isdisjoint(target_ids):
        mismatch_ids_warning = "There are no " + " or ".join(feature_types) + " features with matching IDs in the " \
                                                                         "reference " \
                                                                     "and target annotation"
        warnings.warn(mismatch_ids_warning)


if __name__ == "__main__":
    main()


