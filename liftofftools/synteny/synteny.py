import numpy as np
from liftofftools.cli_arguments import ARGS
from liftofftools.synteny import plot_gene_order
from liftofftools import filepaths, annotation
import nltk




def analyze_synteny(ref_db, target_db,  ref_fa, target_fa):
    print('Analyzing synteny')
    output_file = filepaths.build_filepath([ARGS.dir, filepaths.SYNTENY_OUTPUTS['gene_order']])
    if filepaths.make_file(output_file):
        less_contigs_fa, more_contigs_fa, less_contigs_db, more_contigs_db = find_most_contiguous_genome(ref_fa,
                                                                                                      target_fa,
                                                                                              ref_db, target_db)
        default_chrom_order = order_chroms_by_size(less_contigs_fa)
        default_gene_order = order_genes(default_chrom_order, less_contigs_db)
        matched_chrom_order = match_chrom_order(default_gene_order, more_contigs_db)
        matched_gene_order = order_genes(matched_chrom_order ,more_contigs_db)
        ref_gene_order, target_gene_order= get_ref_and_target_gene_order(default_gene_order, matched_gene_order,
                                                                         ref_fa==more_contigs_fa)
        output_rows = print_order_output(ref_gene_order, target_gene_order, ref_db,target_db, output_file)
        plot_gene_order.plot_gene_order(output_rows)

def find_most_contiguous_genome(ref_fa, target_fa, ref_db, target_db):
    if len([seq for seq in ref_fa]) < len([seq for seq in target_fa]):
        return ref_fa, target_fa, ref_db, target_db
    return target_fa, ref_fa, target_db, ref_db


def get_ref_and_target_gene_order(default_gene_order, matched_gene_order, ref_has_more_contigs):
    if ref_has_more_contigs:
        return matched_gene_order, default_gene_order
    return default_gene_order, matched_gene_order



def order_chroms_by_size(fa):
    chroms_with_size = []
    for chrom, seq in fa.items():
        chroms_with_size.append([chrom, len(seq)])
    chroms_with_size.sort(key=lambda x: x[1], reverse=True)
    ordered_chroms = np.array(chroms_with_size)[:,0]
    return get_chrom_order_dict(ordered_chroms)


def get_chrom_order_dict(ordered_chroms):
    chrom_order_dict = {}
    for i in range(len(ordered_chroms)):
        chrom_order_dict[ordered_chroms[i]] = i
    return chrom_order_dict



def match_chrom_order(ref_gene_order, target_db):
    chrom_to_ref_gene_order = get_default_order_of_target_genes(target_db, ref_gene_order)
    target_positions = get_median_default_gene_order(chrom_to_ref_gene_order)
    target_positions.sort(key = lambda x: x[1])
    ordered_chroms = np.array(target_positions)[:,0]
    return get_chrom_order_dict(ordered_chroms)


def get_default_order_of_target_genes(target_db, ref_gene_order):
    chrom_to_ref_gene_order = {}
    all_genes = target_db.get_features_of_type(ARGS.ft)
    for gene in all_genes:
        if gene.seqid not in chrom_to_ref_gene_order:
            chrom_to_ref_gene_order[gene.seqid] = [0]
        if gene.id in ref_gene_order:
            chrom_to_ref_gene_order[gene.seqid].append(ref_gene_order[gene.id])
    return chrom_to_ref_gene_order


def get_median_default_gene_order(chrom_to_ref_gene_order):
    chrom_and_avg_gene_number = []
    for chrom, gene_order_numbers in chrom_to_ref_gene_order.items():
        chrom_and_avg_gene_number.append([chrom, np.median(gene_order_numbers)])
    return chrom_and_avg_gene_number



def order_genes(chrom_order, gene_db):
    gene_order_dict = {}
    all_genes = gene_db.get_features_of_type(ARGS.ft)
    all_genes.sort(key=lambda x: [chrom_order[x.seqid], x.start])
    for i in range(len(all_genes)):
        gene_order_dict[all_genes[i].id] =i
    return gene_order_dict



def print_order_output(ref_gene_order, target_gene_order, ref_db, target_db, output_file):
    if ARGS.edit_distance:
        edit_distance = nltk.edit_distance(list(ref_gene_order.values()), list(target_gene_order.values()),
                                           transpositions=True)
    else:
        edit_distance = 'not calculated'
    output_rows = []
    target_gene_dict = target_db.get_feature_dict(ARGS.ft)
    target_idx_to_output_rows = add_ref_genes_to_print(output_rows, ref_gene_order, target_gene_order, ref_db, target_gene_dict)
    insert_novel_target_genes_to_print(target_idx_to_output_rows, ref_gene_order, target_gene_order,  output_rows,
                                       target_gene_dict)
    with open(output_file, 'w') as f:
        f.write("Edit distance=" + str(edit_distance) + "\n" + "\n")
        for row in output_rows:
            f.write("\t".join([ str(col) for col in row]) + "\n")
    if len(output_rows) == 0:
        sys.exit("No features with matching IDs to plot")
    return output_rows


def add_ref_genes_to_print(output_rows, ref_gene_order, target_gene_order, ref_db, target_gene_dict):
    target_idx_to_output_row = {}
    output_row_idx = -1
    ref_gene_dict = ref_db.get_feature_dict(ARGS.ft)
    for gene in ref_gene_order:
        output_row_idx += 1
        ref_seq_name = ref_gene_dict[gene].seqid
        target_order, target_seq_name = 'NA', 'NA'
        seq_id = -1
        if gene in target_gene_order:
            target_order = target_gene_order[gene]
            target_idx_to_output_row[target_order] =output_row_idx
            target_seq_name = target_gene_dict[gene].seqid
            seq_id = annotation.get_perc_id(target_gene_dict[gene])
        output_rows.append([gene, ref_gene_order[gene], target_order, ref_seq_name, target_seq_name, seq_id])
    return target_idx_to_output_row


def insert_novel_target_genes_to_print(target_idx_to_output_rows, ref_gene_order, target_gene_order, output_rows, target_gene_dict):
    for gene in target_gene_order:
        if gene not in ref_gene_order:
            seq_name = target_gene_dict[gene].seqid
            seq_id = annotation.get_perc_id(target_gene_dict[gene])
            if target_gene_order[gene] == 0:
                previous_ind = -1
            else:
                previous_ind = target_idx_to_output_rows[target_gene_order[gene]-1]
            output_rows.insert(previous_ind +1, [gene, 'NA', target_gene_order[gene], 'NA', seq_name, seq_id])
            target_idx_to_output_rows[target_gene_order[gene]] = previous_ind + 1


