To run Epsilon, use the following command:

```bash
epsilon --config <name of config file> --jobs <number of parallel jobs>
```

The ```config``` file defines all the input and output files as well as all the parameters of Epsilon.
The structure is the following:

```bash
reference_genome:
negative_control_bamlist:
loci_list:
alt_mode: "generic"
nranks: 1
max_noise_level: 0.05
tumor_path:
matched_normal_path:
calling_loci_list:
variant_calling_model: "binomial"
fdr_method: "BH"
alpha: 0.05
prior: 0.3
posterior_cutoff: 0.5
output_vcf: "output.vcf"
```

For a detailed description of all the elements of the config file please look at the ```Config file structure``` subpage.