import pandas as pd
import pyBigWig
import numpy as np
import math

atac_path1 = 'hg19_ATAC/hg19_t0_atac.bw'
atac_path2 = 'hg19_ATAC/hg19_t24h_atac.bw'

snp_position = 61595564

bw = pyBigWig.open(atac_path1)
max_range = 10_000
data = bw.values("chr11", snp_position - max_range - 1, snp_position + max_range)
data = [0 if math.isnan(x) else x for x in data]
data = np.array(data)
np.save('hg19_ATAC/rs968567_small_h0_atac.npy',data)
max_range = 100_000
data = bw.values("chr11", snp_position - max_range - 1, snp_position + max_range)
data = [0 if math.isnan(x) else x for x in data]
data = np.array(data)
np.save('hg19_ATAC/rs968567_large_h0_atac.npy',data)

bw = pyBigWig.open(atac_path2)
max_range = 10_000
data = bw.values("chr11", snp_position - max_range - 1, snp_position + max_range)
data = [0 if math.isnan(x) else x for x in data]
data = np.array(data)
np.save('hg19_ATAC/rs968567_small_h24_atac.npy',data)
max_range = 100_000
data = bw.values("chr11", snp_position - max_range - 1, snp_position + max_range)
data = [0 if math.isnan(x) else x for x in data]
data = np.array(data)
np.save('hg19_ATAC/rs968567_large_h24_atac.npy',data)
