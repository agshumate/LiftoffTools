from liftofftools import alignment



class Cluster():
    def __init__(self, rep, members=[]):
        self.cluster_rep = rep
        self.members = members

    def add_member(self, member):
        self.members.append(member)

    def remove_members(self, members):
        for member in members:
            self.members.remove(member)

    def __str__(self):
        return ",".join(self.members)



def make_clusters_dict(tsv_file):
    all_clusters = {}
    if tsv_file is not None:
        with open(tsv_file) as tsv:
            for line in tsv:
                rep, member = line.strip().split("\t")
                if rep in all_clusters:
                    all_clusters[rep].add_member(member)
                else:
                    all_clusters[rep] = Cluster(rep, [member])
    return all_clusters


def get_cluster_member_to_rep(cluster_dict):
    cluster_member_to_rep = {}
    for rep, Cluster in cluster_dict.items():
        for member in Cluster.members:
            cluster_member_to_rep[member] = Cluster.cluster_rep
    return cluster_member_to_rep


def process_unmapped_genes(clusters, gene_db, feature_types):
    unmapped_gene_to_cluster = {}
    gene_ids = gene_db.get_all_parent_feature_ids(feature_types)
    for cluster_rep, cluster in clusters.items():
        members_to_remove = []
        for member in cluster.members:
            if member not in gene_ids:
                unmapped_gene_to_cluster[member] = [gene for gene in cluster.members if gene in gene_ids]
                members_to_remove.append(member)
        cluster.remove_members(members_to_remove)
    return unmapped_gene_to_cluster


def get_seq_ids_of_cluster_mates(gene, cluster_members, ref_seqs, target_seqs):
    seq_ids_of_cluster_mates = {}
    ref_seq = ref_seqs.get_longest_isoform_dict([gene])[gene]
    target_cluster_seq_dict = target_seqs.get_longest_isoform_dict(cluster_members)
    for member, target_seq in  target_cluster_seq_dict.items():
        if ref_seqs.is_protein:
            aln = alignment.ProteinAlignment(ref_seq, target_seq)
        else:
            aln = alignment.DNAAlignment(ref_seq, target_seq)
        seq_ids_of_cluster_mates[member] = float(aln.calculate_seq_id())
    return seq_ids_of_cluster_mates



def write_cluster_output(cluster_group1, cluster_group2, output_file, is_coding):
    cluster_id = 0
    if is_coding:
        prefix = 'protein'
    else:
        prefix = 'DNA'
    for rep, cluster1 in cluster_group1.items():
        cluster_id += 1
        cluster_label = prefix +"_cluster_" + str(cluster_id)
        size_cluster_1 = len(cluster1.members)
        cluster2 = cluster_group2[rep]
        size_cluster2 = len(cluster2.members)
        output_file.write("\t".join([cluster_label, str(size_cluster_1), str(cluster1), str(size_cluster2), str(cluster2)]) +
                "\n")