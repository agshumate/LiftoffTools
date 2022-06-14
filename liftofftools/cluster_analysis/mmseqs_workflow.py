from liftofftools.subprocesses import SubprocessCaller
from liftofftools.cli_arguments import ARGS
from liftofftools.filepaths import make_file, remove_file


def new_cluster_workflow(fasta_file_list, output_prefix):
    mmseqs_path = get_mmseqs_path()
    db_name = create_db(mmseqs_path, fasta_file_list, output_prefix)
    cluster_name = cluster(mmseqs_path, db_name)
    tsv_name = create_tsv(mmseqs_path, db_name, cluster_name)
    return tsv_name


def get_mmseqs_path():
    if ARGS.mmseqs_path is None:
        return 'mmseqs'
    return ARGS.mmseqs_path

def create_db(mmseqs_path, fasta_file_list, output_name):
    subcommand = 'createdb'
    positionals = fasta_file_list + [output_name + "_db"]
    command_line_args = {'-v': '2'}
    if make_file(output_name + "_db"):
        SubprocessCaller(mmseqs_path, subcommand, positionals, command_line_args).run_command()
    return output_name + "_db"



def cluster(mmseqs_path, db_name):
    subcommand = 'linclust'
    positionals = [db_name, db_name+"_cluster", ARGS.dir+ "/mmseqs_intermediates/"+ 'tmp']
    command_line_args = parse_params()
    if make_file(db_name + "_cluster.dbtype"):
        remove_file(db_name + "_cluster.dbtype")
        SubprocessCaller(mmseqs_path, subcommand, positionals, command_line_args).run_command()
    return db_name+"_cluster"



def create_tsv(mmseqs_path, db_name, cluster_name):
    subcommand = 'createtsv'
    positionals = [db_name, db_name, cluster_name, cluster_name + ".tsv"]
    command_line_args = {'-v': '2'}
    if make_file(cluster_name + ".tsv"):
        SubprocessCaller(mmseqs_path, subcommand, positionals, command_line_args).run_command()
    return cluster_name + ".tsv"



def parse_params():
    command_line_args = {'-v': '2'}
    split_args = ARGS.mmseqs_params.split(" ")
    for i in range(0, len(split_args)-1,2):
        command_line_args[split_args[i]] = split_args[i+1]
    if "-c" not in command_line_args:
        command_line_args['-c'] = '0.9'
    if '--min-seq-id' not in command_line_args:
        command_line_args['--min-seq-id']='0.9'
    return command_line_args
