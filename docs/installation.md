
# Installation

Install Epsilon from its repository:

```bash
git clone https://github.com/JakubLiu/Epsilon.git
```

After cloning the repository the following directory structure is present:

```text
├── env.yaml
└── project/
    ├── Epsilon.smk
    ├── Epsilon_config.yaml
    ├── run_epsilon.sh
    └── scripts/
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

In order to satisfy all dependencies please execute the following commands:

```bash
cd Epsilon
conda env create -f env.yaml
conda activate Epsilon_env
```