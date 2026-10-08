
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