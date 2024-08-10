# Methven
Predicting the effect of non-coding mutations on single-cell DNA methylation using DNA language model.

<p align="center"><img width="100%" src="model/model.png" /></p>

## Environment
- Python == 3.9
- Tensorflow-gpu == 2.7
- Protobuf == 3.20
- Scikit-learn == 1.1

## Training
You can specify the model size and other hyper-parameters through the command:
```shell
cd [work_path]
python training.py -m small --epoch 100 --lr 0.005 --save_dir model/weights/
```
Note: 
- When using the preprocessing scripts and the parameters of Methven to process and predict on your own dataset, place all **Reference Allele** to the "before mutation" input branch for hg19 sequence alignment.

## Processed training dataset and model weights
|model|trained weights|parameter|processed training data|test data|
|:---:|:---:|:---:|:---:|:---:|
|Small|[Download](https://www.psymukb.net:83/EMO_Download/trained_weights/small/)|722,786|[Download](https://www.psymukb.net:83/EMO_Download/training_test_set/small/train_small_post.pkl)|[Download](https://www.psymukb.net:83/EMO_Download/training_test_set/small/test_small_post.pkl)|
|Large|[Download](https://www.psymukb.net:83/EMO_Download/trained_weights/large/)|814,946|[Download](https://www.psymukb.net:83/EMO_Download/training_test_set/large/train_large_post.pkl)|[Download](https://www.psymukb.net:83/EMO_Download/training_test_set/large/test_large_post.pkl)|

## Raw data
|Data|resource|
|:---:|:---:|
|single-cell meQTL(CD4+ T cell, Monocyte)|[EPIGEN MeQTL Database](https://epicmeqtl.kcl.ac.uk/ )|
|single-cell ATAC-seq(CD4+ T cell, Monocyte)|[EpiMap](https://personal.broadinstitute.org/cboix/epimap/metadata/Short_Metadata.html)|
|GRCh37/hg19 genome|[UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway)|
|Regulation annotation (hg19)|[UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19)|
|Unstimulated and 24hr stimulated ATAC-seq|[GWAS summary statistics](http://plaza.umin.ac.jp/~yokada/datasource/software.htm)|
|RA-associated SNPs|[GEO Series accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138767)|