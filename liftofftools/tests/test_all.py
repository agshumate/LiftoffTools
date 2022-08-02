import pytest
import liftofftools.liftofftools

def test_yeast():

    #inputs
    ref_fa = 'liftofftools/tests/GCA_000146045.2_R64_genomic.fna'
    target_fa = ref_fa
    ref_gff = 'liftofftools/tests/GCA_000146045.2_R64_genomic.gff'
    target_gff = 'liftofftools/tests/GCA_000146045.2_R64_genomic_out.gff'

    #outputs
    synteny_out = 'liftofftools/tests/liftofftools_output/gene_order'
    variants_out = 'liftofftools/tests/liftofftools_output/variant_effects'
    clusters_out = 'liftofftools/tests/liftofftools_output/clusters'
    unmapped_out = 'liftofftools/tests/liftofftools_output/unmapped_closest_paralogs'

    #expected_outputs
    expected_output_dir = 'liftofftools/tests/liftofftools_output_expected/'
    synteny_exp = expected_output_dir + 'gene_order_expected'
    variants_exp = expected_output_dir + 'variant_effects_expected'
    clusters_exp = expected_output_dir + 'clusters_expected'
    unmapped_exp = expected_output_dir + 'unmapped_closest_paralogs_expected'

    #run
    args = ['all','-r', ref_fa, '-t', target_fa, '-rg', ref_gff, '-tg', target_gff, '-force']
    liftofftools.liftofftools.main(args)


    compare_files(synteny_out, synteny_exp)
    compare_files(variants_out, variants_exp)
    compare_files(clusters_out, clusters_exp)
    compare_files(unmapped_out, unmapped_exp)


def compare_files(observed, expected):
    with open(observed) as fh1, open(expected) as fh2:
        fh1_lines = fh1.readlines()
        fh2_lines = fh2.readlines()
        for i in range(len(fh1_lines)):
            expected_output = fh2_lines[i].strip()
            observed_output = fh1_lines[i].strip()
            assert observed_output == expected_output

