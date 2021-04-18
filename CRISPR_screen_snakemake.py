"""
Analyze CRISPR screening data using bbmerge and guide counting
"""

__pipeline__ = 'CRISPR_screen_snakemake'
__version__ = '0.1'
__author__ = 'Alex Bartlett'


# load standard modules
from glob import glob
from os import path
from datetime import datetime
import argparse
import sys
import re

parent_path = path.dirname(path.dirname(path.realpath(__file__)))
print(parent_path)
sys.path.append(parent_path)
from src.utils import pipelining, registry, s3

try:
    import pandas as pd
    from snakemake import snakemake
    from multiprocessing import cpu_count
    import yaml
    import psutil

except ModuleNotFoundError as e:
    print(e)
    print('===========')
    print('Failed to find a required module. Did you forget to activate the "ALB" environment?')
    print('Try running:\n conda activate ALB')
    exit()

PLASMID_LIB = '200424_plasmid_norms_RCed.csv'

def parse_args():
    parser = argparse.ArgumentParser(description="CRISPR screening analysis pipeline",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = pipelining.add_io_args(parser, default_strand="stranded")
    parser = pipelining.add_qc_args(parser)
    parser = pipelining.add_snakemake_args(parser)
    parser = pipelining.add_arrakis_args(parser)
    parser.add_argument('--mutant_table',   help='FASTA file containing mutant sequences')
    parser.add_argument('--offset',       help='Offset from the start of the read which contains the mutagenized sequence. Unused if --naive is set.', type=int, default=23)
    parser.add_argument('--naive',        help='Brute force lookup of all sequences in each read. Slow.', action='store_true')
    parser.add_argument('--control_sgrnas', type = str, help='List of sgRNAs to be used as negative control set for MaGECK null distribution')
    parser.add_argument('--mageck_norm_method', help="which group to normalize counts to. default is entire library but if 'control' is used, only negative control sgRNAs will be included")
    parser.add_argument('--vs_plasmid_lib', nargs = '?', help="analyze vs. plasmid library guide abundances instead of unsorted controls", default=PLASMID_LIB)
    args = parser.parse_args()

    return args


if __name__ == "__main__":

    # ------------------------------------
    # Prepare inputs and parameters
    # ------------------------------------
    args = parse_args()

    # set up config dict variables, set to null if not present, required for the snakefile to not complain
    config_dict = {'version': __version__,
                   'SAMPLE_SHEET': path.realpath(args.sample_sheet) if args.sample_sheet else "",
                   "CORES": args.num_cores,
                   "MUTANT_TABLE": args.mutant_table,
                   "OFFSET": args.offset,
                   "NAIVE": args.naive,
                   "CONTROL_SGRNAS": args.control_sgrnas, 
                   "MAGECK_NORM_METHOD": args.mageck_norm_method,
                   "VS_PLASMID_LIB": path.realpath(args.vs_plasmid_lib) if args.vs_plasmid_lib else ""
                   }

    # parse the ASR here to use for outprefix if necessary
    parsed_asr = registry.parse_ASR_number(path.realpath(args.output_dir))
    if args.asr:
        assert re.match(r'^ASR\-\d{4}$', args.asr), "Please specify asr as ASR-[4 digits], e.g. ASR-0010"
        if parsed_asr:
            assert args.asr == parsed_asr, f"Error: Specified {args.asr}, but parsed {parsed_asr}"
            asr = args.asr
        else:
            asr = args.asr
    else:
        asr = parsed_asr

    if args.outprefix:
        config_dict["OUTPREFIX"] = args.outprefix
        outprefix = args.outprefix
    else:
        config_dict["OUTPREFIX"] = asr
        outprefix = asr

    # link input reads and store in config dict
    fastq_dict = pipelining.link_input_reads(args)
    config_dict["NAMES"] = [x for x in fastq_dict.keys()]

    # extract sample groups from sample sheet
    sample_sheet_df = pipelining.load_dataframe(config_dict["SAMPLE_SHEET"])
    cells = sample_sheet_df.Group.unique()
    unique_groups = list()
    for cell in cells:
        groups = cell.split(', ')
        for group in groups:
            if not(group in unique_groups):
                unique_groups.append(group)

    config_dict["GROUPS"] = unique_groups

    # extract number of control samples and experimental samples per group and use it to generate mageck args
    for group in unique_groups:
        controls = sample_sheet_df[sample_sheet_df['Group'].str.contains(group) & sample_sheet_df['Control'].str.contains('Yes')].Name
        samples = sample_sheet_df[sample_sheet_df['Group'].str.contains(group) & sample_sheet_df['Control'].str.contains('No')].Name
        num_controls = len(controls)
        num_samples = len(samples)

        control_arg = ""
        sample_arg = ""

        if (config_dict['VS_PLASMID_LIB'] is not None):
            assert num_controls == 0, f"Error: Controls were specified in sample sheet but --vs_plasmid_lib was also included"
            control_arg = 0
            for i in range(1, num_samples):
                sample_arg = sample_arg + str(i) + ','
            sample_arg = sample_arg[:-1]
            
        else:
            assert num_controls != 0, f"Error: No controls specified in sample sheet but --vs_plasmid_lib not included"
            for i in range(0, num_controls + num_samples):
                if (i < num_controls):
                    control_arg = control_arg + str(i) + ','
                else:
                    sample_arg = sample_arg + str(i) + ','
            control_arg = control_arg[:-1]
            sample_arg = sample_arg[:-1]

        config_dict["CONTROL_ARG_" + str(group)] = control_arg
        config_dict["SAMPLE_ARG_" + str(group)] = sample_arg

    # set up numbering and naming convention for results file structure
    prestitch_dir = "0_cleaned_raw"
    prestitch_dir_num = 0

    config_dict["PRESTITCH_DIR"] = prestitch_dir
    config_dict["PRESTITCH_DIR_NUM"] = prestitch_dir_num

    premap_dir = f"{prestitch_dir_num+1}_stitched"
    config_dict["PREMAP_DIR"] = premap_dir

    premap_dir_num = prestitch_dir_num + 1

    config_dict["PREMAP_DIR_NUM"] = premap_dir_num
    postmap_dir_num = premap_dir_num + 1
    config_dict["POSTMAP_DIR_NUM"] = postmap_dir_num
    postmap_dir = f"{postmap_dir_num}_counts"
    config_dict["POSTMAP_DIR"] = postmap_dir

    # set string to be used as end of filename of count tables
    if config_dict["NAIVE"]:
        name_end_str="naive"
    else:
        name_end_str=f"offset_{args.offset}"

    config_dict["NAME_END_STR"] = name_end_str

    # dump config file to yaml
    if args.config_file:
        config_file = args.config_file
    else:
        config_file = path.join(args.output_dir, outprefix+'_config.yaml')
    config_FH = open(config_file, 'w')
    yaml.dump(config_dict, config_FH, default_flow_style=False)
    config_FH.close()

    # Define final output (snakemake targets)
    targets = []

    for name in fastq_dict:

        # QC outputs
        if args.fastqc:
            targets.append(f"1_fastqc/{name}_R1_fastqc.html")
            targets.append(f"1_fastqc/{name}_R2_fastqc.html")

        # stitched reads and count table outputs
        targets.append(f"{premap_dir}/{name}.stitched.fastq")
        targets.append(f"{postmap_dir}/{name}.{name_end_str}.counts.txt")

    # mageck input and output files for each group of samples to be averaged and/or normalized together
    for group in unique_groups:
        targets.append(f"3_mageck_inputs/{group}.txt")

        if ((config_dict["MAGECK_NORM_METHOD"] is not None) & (config_dict["CONTROL_SGRNAS"] is not None)):
            targets.append(f"3_mageck/{group}/{group}_vs_control_sgRNAs.gene_summary.txt")
        else:
            targets.append(f"3_mageck/{group}/{group}.gene_summary.txt")


    # with the outprefix defined, save the command ran
    pipelining.save_snakemake_cmd(args.output_dir, outprefix, sys.argv)

    if not args.create_dag and args.print_targets:
        print(targets)

    # # Check resource availability
    # mem = psutil.virtual_memory().available / 1000000

    # Run pipeline with snakemake
    snakemake(path.join(path.dirname(path.dirname(path.realpath(__file__))), "snakemake", "CRISPR_screen.smk"),
              targets=targets, forcetargets=args.overwrite, keepgoing=args.keepgoing,
              workdir=args.output_dir, cores=args.num_cores, dryrun=args.dryrun, config=config_dict,
              restart_times=args.restart_times, force_incomplete=args.force_incomplete, unlock=args.unlock,
              notemp=args.keep_temp, printdag=args.create_dag)

if not args.create_dag:
    print('======================')

if not args.dryrun and not args.create_dag and args.asr:

    analysis_type = "CRISPR screening analysis pipeline"
    if args.registry:
        bucket_link = registry.get_s3_bucket(asr)
        print(bucket_link)
        registry.track_seqRunAnalysis(asr, analysis_type, __pipeline__, __version__, sys.argv, bucket_link)
    if args.s3:
        bucket_link = registry.get_s3_bucket(asr)
        s3.sync_to_s3(args.output_dir, bucket_link)

if not args.create_dag:
    print('======================')

