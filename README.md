# Methven
Predicting the effect of non-coding mutations on single-cell DNA methylation using DNA language model.

<p align="center"><img width="100%" src="model/model.png" /></p>

## Environment
- Python == 3.9
- Tensorflow-gpu == 2.7
- Protobuf == 3.20
- Scikit-learn == 1.1

## Running interface
Please download the reference genome and the pretrained model weights from the [Cloud Storage](https://www.psymukb.net:83/EMO_Download/trained_weights/) (download all weights and save to one folder). Please keep the same file name as when you downloaded it, and the program will automatically identify which model to use. 

```python
import numpy as np
from src.utils_sign_prediction import *
from src.utils_slope_prediction import *


```

## Training
You can specify the model size and other hyper-parameters through the command:
```shell
cd [work_path]
python training.py -m small --epoch 100 --lr 0.005 --save_dir model/weights/
```
Note: 
- When using the preprocessing scripts and the parameters of Methven to process and predict on your own dataset, place all **Reference Allele** to the "before mutation" input branch for hg19 sequence alignment.

## Processed training dataset and model weights
|model|trained weights|slope trained weights|parameter|processed training data|test data|
|:---:|:---:|:---:|:---:|:---:|
|Small|[Download](https://www.psymukb.net:83/Methven_Download/small/weight/)|722,786|[Download](https://www.psymukb.net:83/Methven_Download/small/weight_slope/)|[Download](https://www.psymukb.net:83/Methven_Download/small/data/small_train.dataset)|[Download](https://www.psymukb.net:83/Methven_Download/small/data/small_test.dataset)|
|Large|[Download](https://www.psymukb.net:83/Methven_Download/large/weight/)|814,946|[Download](https://www.psymukb.net:83/Methven_Download/large/weight_slope/)|[Download](https://www.psymukb.net:83/Methven_Download/large/data/large_train.dataset)|[Download](https://www.psymukb.net:83/Methven_Download/large/data/large.dataset)|

## Raw data
|Data|resource|
|:---:|:---:|
|single-cell meQTL(CD4+ T cell, Monocyte)|[EPIGEN MeQTL Database](https://epicmeqtl.kcl.ac.uk/ )|
|single-cell ATAC-seq(CD4+ T cell, Monocyte)|[EpiMap](https://personal.broadinstitute.org/cboix/epimap/metadata/Short_Metadata.html)|
|GRCh37/hg19 genome|[UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway)|
|Regulation annotation (hg19)|[UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19)|
|Unstimulated and 24hr stimulated RA ATAC-seq|[GWAS summary statistics](http://plaza.umin.ac.jp/~yokada/datasource/software.htm)|
|RA-associated SNPs|[GEO Series accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138767)|
