Epsilon can be used to call somatic point mutations directly from the circulating tumor DNA.

The steps are as follows:

1. Fit and error model on a negative control cohort (ideally this cohort would comprise of normal samples and not liquid biopsy samples)
2. Use the fitted error model to call SNV's in the cell-free DNA (ctDNA) sample

**important**
The variant calling should not be done on the consensus BAM file, but on the pre-consensus BAM. For more information on the pipeline used
to create the alignment files please reffer to: ```https://github.com/fulcrumgenomics/fgbio```.

The input data and parameters are specified in a ```config``` file of the structure:

```bash
# Configuration file for running Epsilon

# path to the reference genome
reference_genome:

# path to the text file containing paths to the bam files that make up the negative control cohort
negative_control_bamlist:

# path to the csv file containing the list of loci on which the error model will be fitted
loci_list:

# mode: 'generic' or 'specific' (default 'generic')
alt_mode: "generic"

# number of ranks (default 1)
nranks: 1

# maximum noise level, above which a locus is excluded (default 0.05)
max_noise_level: 0.05

# path to the text file that contains the path to the tumor bam file
tumor_path:


# path to the csv file that defines the loci in which we want to look for SNVs
calling_loci_list:

# variant calling mode 'bayesian' or 'binomial' (default 'binomial')
variant_calling_model: "binomial"

# FDR method, only applicable if the variant_calling_model is set to 'binomial'. Possible values: ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"] (default 'BH')
fdr_method: "BH"

# signficance level, only applicable if the variant_calling_model is set to 'binomial' (default 0.05)
alpha: 0.05

# Bayesian prior, only applicable if the variant_calling_model is set to 'bayesian' (default 0.3)
prior: 0.3

# Bayesian posterior cutoff, only applicable if the variant_calling_model is set to 'bayesian' (default 0.5)
posterior_cutoff: 0.5

# name of the output vcf file
output_vcf: "output.vcf"
```

For more information please look at the ```Usage/'Config file structure'``` subpage.

Next, to call variants, please execute:

```bash
epsilon_cfDNA --config <config file> --jobs <number of parallel jobs>
```

**important**

cfDNA data is usually panel data, i.e. the coverage is restricted to a relatively small region. Therefore the user should
be careful when defining the ```calling_loci_list``` argument in the ```config``` file, so that it is contained within the panel regions.
If no information about the panel regions is available, a workaround would be to use ```samtools depth``` to infer the panel. (```https://www.htslib.org/doc/samtools-depth.html```).