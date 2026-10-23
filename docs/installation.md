
# Installation

Install Epsilon from its repository:

```bash
git clone https://github.com/JakubLiu/Epsilon.git
```

In order to satisfy all dependencies please execute the following commands:

```bash
cd Epsilon
conda env create -f env.yaml
conda activate Epsilon_env
```

In addition add the executable to the ```PATH``` variable (assuming you are located in the ```Epsilon/``` directory):

```bash
cd code
export PATH="$PWD:$PATH"
```
Optionally you can make it permanent by applying the change to your ```.bashrc``` file:

```bash
echo 'export PATH="/absolute/path/to/your/directory:$PATH"' >> ~/.bashrc
```

Where ```/absolute/path/to/your/directory``` is the absolute path to ```Epsilon/code```.

The final step is to grant yousrself execution rights to all scripts. Assuming you are located in ```Epsilon/code/``` run the following:

```bash
chmod +x *
chmod +x scripts/*
```

To test the installation type:

```bash
epsilon --help
```

The following output should be printed to the standard output:

```bash
Usage: epsilon [OPTIONS]

Run the Epsilon variant-calling pipeline using Snakemake.

Options:
  --config FILE   Path to the Epsilon configuration YAML file
  --jobs N        Number of CPU cores to use
  -h, --help      Show this help message and exit

Example:
  epsilon --config Epsilon_config.yaml --jobs 8
```
