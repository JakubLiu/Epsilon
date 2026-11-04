You do not always have to run the entire Epsilon workflow if you already have a fitted error model and just want to
apply it to a new sample (or the same sample at different loci).
In that case the ```epsilon_call``` command comes in handy. The inputs and parameters need to be specified in a ```config``` file
that has the following structure:

```text
reference_genome:
error_model:
alt_mode:
nranks: 1
tumor_path:
matched_normal_path:
calling_loci_list:
variant_calling_model: "binomial"
fdr_method: "BH"
alpha: 0.05
prior: 0.3
posterior_cutoff: 0.5
output_vcf: "output_call.vcf"
```

The only new parameter is the ```error_model``` parameter, that is the path to the fitted ```.rds``` model.
Assuming that you have run the full Epsilon workflow in the ```bioinfo/``` directory, then the error model should be located at
```bioinfo/fitted_error_model/error_model.rds```. Look into the ```"Config file structure"``` subsubpage inside the ```"Usage"``` subpage to a detailed
description of all the parameters of the ```config``` file.

<br>
To run Epsilon call simply execute:


```bash
epsilon_call --config <name of config file> --jobs <number of parallel jobs>
```

### warning
It is important that the ```alt_mode``` used during fitting of the model is the same as the ```alt_mode``` specified in the ```config``` file
that is used to run ```epsilon_call```.