#TENxMultiomeTools

This is the repository of the TENxMultiomeTools an R package which collects multiple tools for Quality Control of scRNAseq and scATACseq data.

Actual installation and loading

```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("drighelli/TENxMultiomeTools")
library(TENxMultiomeTools)
```

