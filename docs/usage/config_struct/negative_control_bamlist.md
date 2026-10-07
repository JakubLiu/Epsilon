### dtype
```str```

### default value
```no default value```

### meaning
Path to the text file that contains the paths to the bam files that are part of the negative control cohort.

### example
Entry in the ```config``` file:
```bash
negative_control_bamlist: home/files/negative_control_list.txt
```

```negative_control_list.txt```:
```text
/home/files/patientA.bam
/home/file/patientsB.bam
/home/files/patientC.bam
/home/file/patientsK.bam
```