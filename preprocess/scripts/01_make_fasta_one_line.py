# make mutation single
import pandas as pd
from pandas import read_parquet
pd.set_option('display.max_columns', None)

fasta_path = '../../datasets/raw/reference_genome_hg19/'

for chr_num in range(1, 23):
    with open(fasta_path + 'chr' + str(chr_num) + '.fa') as fa:
        line = fa.readline()
        line = fa.readline()
        with open(fasta_path + 'chr' + str(chr_num) + '.fasta', 'w+') as fasta:
            while line:
                fasta.write(line.strip())
                line = fa.readline()