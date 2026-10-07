### dtype
```str```

### default value
```no default value```

### meaning
Path to the csv file that defines the loci in which we want to look for SNVs.

### example

Entry in the ```config``` file:
```bash
calling_loci_list: home/files/calling_loci_list.csv
```

Where ```calling_loci_list.csv``` looks like this in the ```"generic"``` ```alt_mode```:

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

### other considerations
Note that the chromosome identifiers must match the ones used in the reference genome fasta file, to which your reads have been aligned.
For example, if the reference uses a ```chr1``` identifier for the 1st chromosome, then the same format needs to be used in your csv file.