##-------------------------------------------------------------
## master snakefile for running CRISPR screen analysis pipeline
##-------------------------------------------------------------

include: "seq_stitching.smk"
include: "seq_lib_mapping.smk"


# generate text files required as inputs by MAGeCK
# for each group of experimental/control samples to be averaged/considered together
rule make_mageck_inputs:
    input:
        expand("{postmap_dir}/{name}.{name_end}.counts.txt", postmap_dir=config["POSTMAP_DIR"],
        name=config["NAMES"], name_end= config["NAME_END_STR"], offset=config["OFFSET"])
    output:
        expand("3_mageck_inputs/{group}.txt", group = config["GROUPS"])
    params:
        sample_sheet = config["SAMPLE_SHEET"],
        mutant_table = config["MUTANT_TABLE"],
        groups = config["GROUPS"],
        vs_plasmid_lib = config["VS_PLASMID_LIB"]
    script:
        "../src/CRISPR_screening/make_mageck_inputs.R"

## run MAGeCK using counts from entire library as null distribution
if ((config["MAGECK_NORM_METHOD"] is not None) & (config["CONTROL_SGRNAS"] is not None)):
	rule run_mageck:
	    input:
	        count_file = "3_mageck_inputs/{group}.txt"
	    output:
	        output_dir = "3_mageck/{group}/{group}_vs_control_sgRNAs.gene_summary.txt"
	    params:
	        control_cols = lambda w: config["CONTROL_ARG_" + w.group],
	        treatment_cols = lambda w: config["SAMPLE_ARG_" + w.group],
	        control_sgRNAs = config["CONTROL_SGRNAS"],
	        mageck_norm_method = config["MAGECK_NORM_METHOD"],
	        out_prefix = "3_mageck/{group}/{group}_vs_control_sgRNAs"

	    shell:
	        "mageck test -k {input.count_file} -c {params.control_cols} -t {params.treatment_cols} --control-sgrna {params.control_sgRNAs} --norm-method {params.mageck_norm_method} -n {params.out_prefix}"
            
## run MAGeCK using counts from only negative control sgRNAs in provided file asa null distribution            
else:
	rule run_mageck:
	    input:
	        count_file = "3_mageck_inputs/{group}.txt"
	    output:
	        output_dir = "3_mageck/{group}/{group}.gene_summary.txt"
	    params:
	        control_cols = lambda w: config["CONTROL_ARG_" + w.group],
	        treatment_cols = lambda w: config["SAMPLE_ARG_" + w.group],
	        out_prefix = "3_mageck/{group}/{group}"

	    shell:
	        "mageck test -k {input.count_file} -c {params.control_cols} -t {params.treatment_cols} -n {params.out_prefix}"