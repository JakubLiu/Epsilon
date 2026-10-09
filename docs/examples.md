### Example 1: full Epsilon workflow

At the start of our analysis, the data is organized in this structure:

```bash
.
├── calling_loci_list.csv
├── loci_list.csv
├── matched_normal_path.txt
├── matched_normal_sample
│   ├── patient_10.sorted.bam
│   └── patient_10.sorted.bam.bai
├── negative_control_bamlist.txt
├── negative_control_samples
│   ├── patient_0.sorted.bam
│   ├── patient_0.sorted.bam.bai
│   ├── patient_1.sorted.bam
│   ├── patient_1.sorted.bam.bai
│   ├── patient_2.sorted.bam
│   ├── patient_2.sorted.bam.bai
│   ├── patient_3.sorted.bam
│   ├── patient_3.sorted.bam.bai
│   ├── patient_4.sorted.bam
│   ├── patient_4.sorted.bam.bai
│   ├── patient_5.sorted.bam
│   ├── patient_5.sorted.bam.bai
│   ├── patient_6.sorted.bam
│   ├── patient_6.sorted.bam.bai
│   ├── patient_7.sorted.bam
│   ├── patient_7.sorted.bam.bai
│   ├── patient_8.sorted.bam
│   ├── patient_8.sorted.bam.bai
│   ├── patient_9.sorted.bam
│   └── patient_9.sorted.bam.bai
├── reference_genome
│   ├── genome.fa
│   ├── genome.fa.amb
│   ├── genome.fa.ann
│   ├── genome.fa.bwt
│   ├── genome.fa.fai
│   ├── genome.fa.pac
│   └── genome.fa.sa
├── tumor_path.txt
└── tumor_sample
    ├── patient_11.bam
    └── patient_11.bam.bai
```

We want to call variants in ```tumor_sample/patient_11.bam```. The matched normal sample is
```matched_normal_sample/patient_10.sorted.bam``` and the samples that form the negative control
cohort are located in the ```negative_control_samples/``` directory. The small reference genome can be found at ```reference_genome/```. All the example data can be downloaded from the GitHub repository.

##### step 1: define the ```config``` file

For this analysis and example of a valid ```config``` file might look like this (assuming that we use the ```bayesian```) variant calling mode:

```text
# Configuration file for running Epsilon

# path to the reference genome
reference_genome: "genomics/reference_genome/genome.fa"

# path to the text file containing paths to the bam files that make up the negative control cohort
negative_control_bamlist: "genomics/negative_control_bamlist.txt"

# path to the csv file containing the list of loci on which the error model will be fitted
loci_list: "genomics/loci_list.csv"

# mode: 'generic' or 'specific' (default 'generic')
alt_mode: "generic"

# number of ranks (default 1)
nranks: 4

# maximum noise level, above which a locus is excluded (default 0.05)
max_noise_level: 0.05

# path to the text file that contains the path to the tumor bam file
tumor_path: "genomics/tumor_path.txt"

# path to the text file that contains the path to the matched normal bam file
matched_normal_path: "genomics/matched_normal_path.txt"

# path to the csv file that defines the loci in which we want to look for SNVs
calling_loci_list: "genomics/calling_loci_list.csv"

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

All the ```.txt``` files are text files with paths to the corresponding ```bam``` files, one path per line.
The ```.csv``` follow the structure as discussed in ```Usage/'Config file structure'/loci_list```.