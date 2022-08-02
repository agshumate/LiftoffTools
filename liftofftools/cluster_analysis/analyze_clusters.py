from liftofftools.cluster_analysis import clusters,mmseqs_workflow
import operator
from liftofftools.filepaths import MMSEQS_INTERMEDIATES as MI, CLSUTER_OUTPUTS, build_filepath, make_directory,\
    remove_directory, make_file, remove_file
import copy


def main(ref_proteins, target_proteins, ref_trans, target_trans, ref_db, target_db, args):
    if args.force:
        remove_directory(args.dir+ "/mmseqs_intermediates")
        remove_file(args.dir + "/clusters")
    mmseqs_dir = build_filepath([args.dir, MI['dir']])
    make_directory(mmseqs_dir)
    run_coding_workflow(ref_proteins, ref_db, target_db, target_proteins,ref_trans, target_trans, args)
    if args.c is False:
        run_noncoding_workflow(ref_db, target_db, ref_trans, target_trans, args)


def run_coding_workflow(ref_proteins, ref_db, target_db, target_proteins,ref_trans, target_trans, args):
    print('Analyzing protein-coding clusters')
    unmapped_coding = analyze_proteins(ref_proteins, ref_db, target_db, target_proteins, args)
    with open(build_filepath([args.dir, CLSUTER_OUTPUTS['unmapped']]), 'w') as f:
        analyze_unmapped(unmapped_coding, ref_trans, target_trans, f, ref_proteins,target_proteins)


def run_noncoding_workflow(ref_db, target_db,ref_trans, target_trans, args):
    print('Analyzing noncoding clusters')
    unmapped_noncoding = analyze_noncoding(ref_trans, ref_db, target_db, target_trans, args)
    with open(build_filepath([args.dir, CLSUTER_OUTPUTS['unmapped']]), 'a') as f:
        analyze_unmapped(unmapped_noncoding, ref_trans, target_trans, f)


def analyze_proteins(ref_proteins, ref_db, target_db, target_proteins, args):
    ref_protein_coding_genes = ref_db.get_protein_coding_features(args.ft)
    target_protein_coding_genes = target_db.get_protein_coding_features(args.ft)
    if len(ref_protein_coding_genes) >0 or len(target_protein_coding_genes) >0:
        ref_clusters, target_clusters = cluster(ref_proteins, target_proteins, ref_protein_coding_genes,
                                                      target_protein_coding_genes,target_db, args)
        with open(build_filepath([args.dir,CLSUTER_OUTPUTS['clusters']]), 'w') as f:
            return process_clusters(ref_clusters, target_clusters, target_db, f, True, args.ft)
    return {}


def cluster(ref_seqs, target_seqs, ref_genes, target_genes, target_db, args):
    ref_output = get_ref_prefix(ref_seqs, args.dir)
    target_output = get_target_prefix(ref_seqs, args.dir)
    ref_tsv = build_clusters(ref_seqs,  ref_genes, ref_output,args)
    target_tsv = build_clusters(target_seqs,  target_genes, target_output, args)
    ref_clusters = clusters.make_clusters_dict(ref_tsv)
    target_clusters = clusters.make_clusters_dict(target_tsv)
    ref_member_to_cluster = clusters.get_cluster_member_to_rep(ref_clusters)
    target_member_to_cluster = clusters.get_cluster_member_to_rep(target_clusters)
    final_target_clusters = add_novel_to_clusters(target_clusters, target_member_to_cluster,
                                                  ref_member_to_cluster, ref_clusters, target_db)
    return ref_clusters, final_target_clusters


def add_novel_to_clusters(target_clusters, target_member_to_cluster, ref_member_to_cluster, ref_cluster, target_db):
    final_target_clusters = copy.deepcopy(ref_cluster)
    for mem in target_member_to_cluster:
        if mem not in ref_member_to_cluster:
            paralog = target_db.get_paralog_name(mem)
            if paralog in ref_member_to_cluster:
                add_to_paralog_cluster(mem, paralog, ref_member_to_cluster, final_target_clusters)
            else:
                cluster_rep = target_member_to_cluster[mem]
                found_cluster = cluster_mates_in_ref_clusters(target_clusters[cluster_rep], ref_member_to_cluster, mem,
                                                              final_target_clusters)
                if found_cluster is False:
                    final_target_clusters[cluster_rep] = target_clusters[cluster_rep]
                    ref_cluster[cluster_rep] = clusters.Cluster(cluster_rep, [])
    return final_target_clusters


def add_to_paralog_cluster(target_gene, paralog, ref_member_to_cluster, final_target_clusters):
    cluster_rep = ref_member_to_cluster[paralog]
    paralog_cluster = final_target_clusters[cluster_rep]
    paralog_cluster.add_member(target_gene)


def cluster_mates_in_ref_clusters(cluster, ref_member_to_cluster, target_gene,final_target_clusters):
    for mate in cluster.members:
        if mate in ref_member_to_cluster:
            cluster_rep = ref_member_to_cluster[mate]
            final_target_clusters[cluster_rep].add_member(target_gene)
            return True
    return False


def analyze_noncoding(ref_trans, ref_db, target_db, target_trans, args):
    ref_noncoding_genes = ref_db.get_noncoding_features(args.ft)
    target_noncoding_genes = target_db.get_noncoding_features(args.ft)
    noncoding_genes = ref_db.get_noncoding_features(args.ft)
    if len(noncoding_genes) >0:
        ref_clusters , target_clusters = cluster(ref_trans, target_trans, ref_noncoding_genes,
                                                       target_noncoding_genes, target_db, args)
        final_cluster_output = build_filepath([args.dir, CLSUTER_OUTPUTS['clusters']])
        with open(final_cluster_output, 'a') as f:
                return process_clusters(ref_clusters, target_clusters, target_db, f, False, args.ft)
    return {}


def build_clusters(sequence_dict, gene_list, ref_output, args):
    fasta_name = ref_output + ".fa"
    if len(gene_list) >0:
        if make_file(fasta_name, args.force):
            longest_seqs = sequence_dict.get_longest_isoform_dict(gene_list)
            write_fasta(longest_seqs, fasta_name)
        tsv = mmseqs_workflow.new_cluster_workflow([fasta_name], ref_output,args)
        return tsv
    return None


def write_fasta(sequence_dict, out_file):
    f = open(out_file, 'w')
    for name, sequence in sequence_dict.items():
        f.write(">" + name + "\n")
        f.write(sequence + "\n")
    f.close


def get_ref_prefix(sequence_dict, dir):
    if sequence_dict.is_protein:
         prefix = MI['ref_coding_prefix']
    else:
        prefix = MI['ref_noncoding_prefix']
    return build_filepath([dir, MI['dir'], prefix])


def get_target_prefix(sequence_dict, dir):
    if sequence_dict.is_protein:
        prefix = MI['target_coding_prefix']
    else:
        prefix = MI['target_noncoding_prefix']
    return build_filepath([dir, MI['dir'], prefix])


def process_clusters(ref_clusters, target_clusters, target_db, output, is_coding, feature_types):
    unmapped_gene_to_cluster = clusters.process_unmapped_genes(target_clusters, target_db, feature_types)
    clusters.write_cluster_output(ref_clusters, target_clusters, output, is_coding)
    return unmapped_gene_to_cluster


def analyze_unmapped(unmapped_genes, ref_trans, target_trans, output_file, ref_proteins={}, target_proteins={}):
    for unmapped, cluster in unmapped_genes.items():

        dna_seq_ids = clusters.get_seq_ids_of_cluster_mates(unmapped, cluster, ref_trans, target_trans)
        if ref_proteins != {}:
            protein_seq_ids = clusters.get_seq_ids_of_cluster_mates(unmapped, cluster, ref_proteins, target_proteins)
            closest_paralog = find_closest_paralog_name(protein_seq_ids)
            protein_seq_id = find_paralog_seqid(closest_paralog, protein_seq_ids)
        else:
            closest_paralog = find_closest_paralog_name(dna_seq_ids)
            protein_seq_id = ''
        dna_seq_id = find_paralog_seqid(closest_paralog, dna_seq_ids)
        output_file.write("\t".join([unmapped, closest_paralog, dna_seq_id, protein_seq_id ]) + "\n")


def find_closest_paralog_name(seq_ids):
    if len(seq_ids) > 0:
        return max (seq_ids.items(), key=operator.itemgetter(1))[0]
    return "None"


def find_paralog_seqid(paralog, seqids):
    if paralog != 'None':
        return str(seqids[paralog])
    else:
        return ''