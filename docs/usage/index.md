In order to run Epsilon, the ```scripts/``` directory must be located in the project directory alongside ```Epsilon.smk```, ```Epsilon_config.yaml``` and 
```run_espilon.sh```. Assuming our project is called ```genome``` a valid directory structure would be:

```bash
└── genome
    ├── Epsilon.smk
    ├── Epsilon_config.yaml
    ├── run_epsilon.sh
    └── scripts
        ├── Epsilon_Fit.sh
        ├── Epsilon_MakeData.sh
        ├── Epsilon_MakeData_call.sh
        ├── Epsilon_call.sh
        ├── call.R
        ├── extract_features_alt_generic.py
        ├── extract_features_alt_generic_call.py
        ├── extract_features_alt_specific.py
        ├── extract_features_alt_specific_call.py
        ├── fit.R
        ├── make_data_alt_generic.py
        ├── make_data_alt_generic_call.py
        ├── make_data_alt_specific.py
        ├── make_data_alt_specific_call.py
        └── remove_germlines.R
```

Assuming that the user is located in ```genome/``` the way to run Epsilon is:

```bash
./run_epsilon.sh --config <name of config file> --jobs <number of parallel jobs>
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