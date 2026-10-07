### detype
```float```

### default value
```0.05```

### meaning
During fitting of the error model, Epsilon tries to remove any mutated sites in order to isolate the error rate.
In an attempt to achieve this, an upper bound is set on the observed variant allele fraction.
The ```max_noise_level``` parameter defines this upperbound.