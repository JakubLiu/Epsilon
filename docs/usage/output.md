After a succesful run, the directory tree will look like this:

```text
project/
├── Epsilon.smk
├── Epsilon_config.yaml
├── run_epsilon.sh
├── scripts/
│   ├── Epsilon_Fit.sh
│   ├── Epsilon_MakeData.sh
│   ├── Epsilon_MakeData_call.sh
│   ├── Epsilon_call.sh
│   ├── call.R
│   ├── extract_features_alt_generic.py
│   ├── extract_features_alt_generic_call.py
│   ├── extract_features_alt_specific.py
│   ├── extract_features_alt_specific_call.py
│   ├── fit.R
│   ├── make_data_alt_generic.py
│   ├── make_data_alt_generic_call.py
│   ├── make_data_alt_specific.py
│   ├── make_data_alt_specific_call.py
│   └── remove_germlines.R
├── data/
├── error_model/
│   └── error_model.rds
├── tumor_data/
│   ├── matched_normal_data.txt
│   └── tumor_data.txt
├── unfiltered_variant_calls/
│   ├── matched_normal.vcf
│   └── tumor.vcf
└── filtered_variant_calls/
    └── tumor.vcf
```

The final output vcf file is the ```filtered_variant_calls/tumor.vcf``` file.