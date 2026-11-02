import os

# set up the paths
EPSILON_DIR = workflow.basedir

SCRIPT_DIR = os.path.join(EPSILON_DIR, "scripts")

MAKE_DATA = os.path.join(SCRIPT_DIR, "Epsilon_MakeData.sh")
MAKE_DATA_CALL = os.path.join(SCRIPT_DIR, "Epsilon_MakeData_call.sh")
FIT_MODEL = os.path.join(SCRIPT_DIR, "Epsilon_Fit.sh")
CALL = os.path.join(SCRIPT_DIR, "Epsilon_call.sh")
REMOVE_GERMLINES = os.path.join(SCRIPT_DIR, "remove_germlines.R")



# get the params from the config
reference_genome = config['reference_genome']
negative_control_bamlist = config['negative_control_bamlist']
loci_list = config['loci_list']
alt_mode = config['alt_mode']
nranks = int(config['nranks'])
max_noise_level = float(config['max_noise_level'])
tumor_sample = config['tumor_path']
matched_normal_sample = config['matched_normal_path']
calling_loci_list = config['calling_loci_list']
variant_calling_model = config['variant_calling_model']
output_vcf = config['output_vcf']


if variant_calling_model == "binomial":
    fdr_method = config['fdr_method']
    alpha = float(config['alpha'])

elif variant_calling_model == "bayesian":
    prior = float(config['prior'])
    posterior_cutoff = float(config['posterior_cutoff'])


rule all:
    input:
        f'filtered_variant_calls/{output_vcf}'



rule make_negative_control_data:
    resources:
        mem_mb = 32000
    input:
        negative_control_bamlist
    output:
        marker = '.snakemake_markers/make_negative_control_data.out',
        data = 'negative_control_data/negative_control_data.txt'
    params:
        ref = reference_genome,
        loci = loci_list,
        bams = negative_control_bamlist,
        mod = alt_mode,
        nr = nranks,
        dirname = 'negative_control_data'
    shell:
        """
        {MAKE_DATA} \
            --bamlist {params.bams} \
            --loci_list {params.loci} \
            --reference_genome {params.ref} \
            --alt_mode {params.mod} \
            --nranks {params.nr} \
            --output_prefix {params.dirname}

        touch {output.marker}
        """


rule fit_error_model:
    resources:
        mem_mb = 32000
    input:
        marker = '.snakemake_markers/make_negative_control_data.out',
        data = 'negative_control_data/negative_control_data.txt'
    output:
        '.snakemake_markers/fit_error_model.out'
    params:
        eps = max_noise_level,
        mod = alt_mode,
        model_name = 'error_model'
    shell:
        """
        {FIT_MODEL} \
            --input {input.data} \
            --noise_level {params.eps} \
            --alt_mode {params.mod} \
            --output_model {params.model_name}

        touch {output}
        """



rule make_tumor_data:
    resources:
        mem_mb = 32000
    input:
        tumor_sample
    output:
        marker = '.snakemake_markers/make_tumor_data.out',
        data = 'tumor_data/tumor_data.txt'
    params:
        ref = reference_genome,
        loci = calling_loci_list,
        mod = alt_mode,
        nr = nranks,
        dirname = 'tumor_data'
    shell:
        """
        {MAKE_DATA_CALL} \
            --bamlist {input} \
            --loci_list {params.loci} \
            --reference_genome {params.ref} \
            --alt_mode {params.mod} \
            --nranks {params.nr} \
            --output_prefix {params.dirname}

        touch {output.marker}
        """

rule make_matched_normal_data:
    resources:
        mem_mb = 32000
    input:
        matched_normal_sample
    output:
        marker = '.snakemake_markers/make_matched_normal_data.out',
        data = 'tumor_data/matched_normal_data.txt'
    params:
        ref = reference_genome,
        loci = calling_loci_list,
        mod = alt_mode,
        nr = nranks,
        dirname = 'matched_normal_data'
    shell:
        """
        {MAKE_DATA_CALL} \
            --bamlist {input} \
            --loci_list {params.loci} \
            --reference_genome {params.ref} \
            --alt_mode {params.mod} \
            --nranks {params.nr} \
            --output_prefix {params.dirname}

        touch {output.marker}
        """

if variant_calling_model == "bayesian":

    rule call_bayesian:
        resources:
            mem_mb = 32000

        input:
            tumor_data = '.snakemake_markers/make_tumor_data.out',
            error_model = '.snakemake_markers/fit_error_model.out'

        output:
            vcf = f'unfiltered_variant_calls/{output_vcf}'

        params:
            prior = float(config["prior"]),
            mod = alt_mode

        shell:
            """
            {CALL} \
                --model fitted_error_model/error_model.rds \
                --input tumor_data/tumor_data.txt \
                --output {output.vcf} \
                --mode bayes_posterior \
                --prior {params.prior} \
                --posterior_cutoff 0.5 \
                --alt_mode {params.mod}
            """

    rule call_bayesian_matched_normal:
        resources:
            mem_mb = 32000

        input:
            matched_normal_data = '.snakemake_markers/make_matched_normal_data.out',
            error_model = '.snakemake_markers/fit_error_model.out'

        output:
            vcf = f'unfiltered_variant_calls/matched_normal_{output_vcf}'

        params:
            prior = float(config["prior"]),
            mod = alt_mode

        shell:
            """
            {CALL} \
                --model fitted_error_model/error_model.rds \
                --input tumor_data/matched_normal_data.txt \
                --output {output.vcf} \
                --mode bayes_posterior \
                --prior {params.prior} \
                --posterior_cutoff 0.5 \
                --alt_mode {params.mod}
            """


elif variant_calling_model == "binomial":

    rule call_binomial:
        resources:
            mem_mb = 32000

        input:
            tumor_data = '.snakemake_markers/make_tumor_data.out',
            error_model = '.snakemake_markers/fit_error_model.out'

        output:
            vcf = f'unfiltered_variant_calls/{output_vcf}'

        params:
            fdr = config["fdr_method"],
            alpha = float(config["alpha"]),
            mod = alt_mode

        shell:
            """
            {CALL} \
                --model fitted_error_model/error_model.rds \
                --input tumor_data/tumor_data.txt \
                --output {output.vcf} \
                --mode binomial_test \
                --fdr_method {params.fdr} \
                --alpha {params.alpha} \
                --alt_mode {params.mod}
            """

    
    rule call_binomial_matched_normal:
        resources:
            mem_mb = 32000

        input:
            matched_normal_data = '.snakemake_markers/make_matched_normal_data.out',
            error_model = '.snakemake_markers/fit_error_model.out'

        output:
            vcf = f'unfiltered_variant_calls/matched_normal_{output_vcf}'

        params:
            fdr = config["fdr_method"],
            alpha = float(config["alpha"]),
            mod = alt_mode

        shell:
            """
            {CALL} \
                --model fitted_error_model/error_model.rds \
                --input tumor_data/matched_normal_data.txt \
                --output {output.vcf} \
                --mode binomial_test \
                --fdr_method {params.fdr} \
                --alpha {params.alpha} \
                --alt_mode {params.mod}
            """


rule remove_germlines:
    resources:
        mem_mb = 32000
    input:
        all_calls = f'unfiltered_variant_calls/{output_vcf}',
        germline_calls = f'unfiltered_variant_calls/matched_normal_{output_vcf}'
    output:
        f'filtered_variant_calls/{output_vcf}'
    shell:
        """
        Rscript {REMOVE_GERMLINES} {input.all_calls} {input.germline_calls} {output}
        """