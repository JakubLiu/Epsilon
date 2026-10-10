### dtype
```str```

### default value
```no default value```

### meaning
Path to the csv file containing the list of loci on which the error model will be fitted.

### example

Entry in the ```config``` file:
```bash
loci_list: home/files/loci_list.csv
```

Where ```loci_list.csv``` looks like this in the ```"generic"``` ```alt_mode```:

```text
"chrom","pos","ref"
"9",133748283,"C"
"9",133738363,"G"
"14",105246551,"C"
"2",29443695,"G"
"2",29432664,"C"
"X",47426121,"C"
"7",140453136,"A"
"7",55259515,"T"
"7",55249071,"C"
```

And like this in the ```"specific"``` ```alt_mode```

```text
"chrom","pos","ref","alt"
"9",133748283,"C","T"
"9",133738363,"G","A"
"14",105246551,"C","T"
"2",29443695,"G","T"
"2",29432664,"C","T"
"X",47426121,"C","G"
"7",140453136,"A","T"
"7",55259515,"T","G"
"7",55249071,"C","T"
```


### Important
The number of loci in the ```loci_list``` file should not be too low. This is due to the fact that the error model needs to be exposed to a sufficient 'sample size', so that the predictor variables exhibit enough variability, and hence parameter estimation is possible. For WES/WGS data and a negative control cohort comprising of 20-100 samples, around 500 loci should be enough. In theory, with a bigger negative control cohort size, fewer loci would be needed, but the cost of including more genomic locations is negligible in terms of compute time. It is therefore advised to maximize their number. This has benefits not only in terms of the parameter identifiability but also the standard error of the parameter estimates.

### other considerations
Note that the chromosome identifiers must match the ones used in the reference genome fasta file, to which your reads have been aligned.
For example, if the reference uses a ```chr1``` identifier for the 1st chromosome, then the same format needs to be used in your csv file.