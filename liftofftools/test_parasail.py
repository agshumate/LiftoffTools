import parasail
import sys
import gffutils

def main():
    ref_fa = sys.argv[1]
    query_fa = sys.argv[2]
    para_ref = parasail.sequences_from_file(ref_fa)
    para_query = parasail.sequences_from_file(query_fa)
    ref_seq = bytes(para_ref[0]).decode('UTF-8').upper()
    query_seq = bytes(para_query[0]).decode('UTF-8').upper()
    gff_file = sys.argv[3]
    matrix_from_filename = parasail.Matrix("planet/matrix")
    feature_db = gffutils.create_db(gff_file, gff_file + "_db", merge_strategy="create_unique",
                                    force=True,
                                    verbose=True)
    # cds_s = feature_db.features_of_type(featuretype="CDS")
    #
    # for cds in cds_s:
    #     ref_seq = ref_seq[:cds.start -1] + ref_seq[cds.start -1: cds.end].replace('A', 'B').replace('C', 'D').replace(
    #         'G', 'E').replace('T', 'F') + ref_seq[cds.end:]
    #
    # print(ref_seq)
    results = parasail.sw_trace(ref_seq, query_seq, 10, 1,
                         matrix_from_filename)
    print(len(results))
    cigar = result.cigar

    print(cigar)
    print(cigar.seq)
    print(cigar.len)
    print(cigar.beg_query)
    print(cigar.beg_ref)
    print(cigar.decode)





main()