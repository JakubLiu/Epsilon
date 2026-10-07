It is possible to submit the different sub-steps of the Epsilon workflow as separate Slurm jobs, in order to speed up computation.

Assuming that a configuration profile is set up as shown in:

```text
https://snakemake.readthedocs.io/en/v7.19.1/executing/cli.html#profiles
```

and assuming that you are located in the ```project/``` directory (and your configuration file is called ```Epsilon_config.yaml```), Epsilon can be submitted as a Slurm job with the following command:

```bash
snakemake --snakefile Epsilon.smk \
          --configfile Epsilon_config.yaml \
          --profile <your profile> \
          --jobs <number of jobs>
```