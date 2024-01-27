import os
import pandas as pd
import numpy as np
import subprocess
from tqdm import tqdm

# step1, clean wrong items, fetch rsids
#'''
file_path = '../../datasets/raw/Kidney_meQTL.txt'
output_path = '../../datasets/middlefile/single_meQTL.txt'
rs_file = '../../datasets/middlefile/meQTL_rsids.txt'
rsid_list = []
with open(file_path) as r:
	with open(output_path,'w+') as w1:
		line = r.readline()
		w1.write(line)
		with open(rs_file,'w+') as w2:
			line = r.readline()
			while line:
				cpg = line.split()[0]
				rsid = line.split()[1]
				A1 = line.split()[2]
				A2 = line.split()[3]
				if(len(A1) == 1 and len(A2) == 1 and cpg != 'NA' and rsid != 'NA'):
					w1.write(line)
					w2.write(rsid + '\n')
				line = r.readline()
#'''
'''
single_meQTL.txt:
CpG	SNP	Ref	Alt	Beta	SE	Pvalue
cg11913416	rs1262461	A	G	-1.179	0.0348	5.6e-119
cg25722041	rs114812780	G	T	1.233	0.0369	3.16e-117
cg25722041	rs12568293	T	C	1.233	0.0369	3.16e-117
cg25722041	rs12117910	A	C	1.233	0.0369	3.17e-117
cg25722041	rs12403339	G	A	1.233	0.0369	3.17e-117
cg25722041	rs12145445	C	G	1.233	0.0369	3.23e-117
cg25722041	rs12120824	C	T	1.233	0.0369	3.24e-117
'''
'''
meQTL_rsids.txt:
rs1262461
rs114812780
rs12568293
rs12117910
rs12403339
'''

# step2 RSID annotation to annovar database
#'''
rs_file = '../../datasets/middlefile/meQTL_rsids.txt'
output_path = '../../datasets/middlefile/dbSNP_rsid_annotation.txt'
vcf_file = '/data/annovar/humandb/hg38_avsnp150.txt'

rs_list = []
with open(rs_file) as r:
	line = r.readline()
	while line:
		rs_list.append(line.strip())
		line = r.readline()
ext_rs_list = list(set(rs_list))
print(len(ext_rs_list)) # 2,810,934

# faster: grep command: grep "\<rs1262461\>" /data/annovar/humandb/hg38_avsnp150.txt -m 1
with open(output_path, 'w+') as w:
	w.write('CHR\tPOS\tRSID\tA1\tA2\n')
	for rsid in ext_rs_list:
		try:
			cmd = 'grep "\<' + rsid + '\>" /data/annovar/humandb/hg38_avsnp150.txt -m 1'
			result = subprocess.check_output(cmd, shell=True) # need decode to string
			SR = result.decode('utf-8')
			CHR = SR.split()[0]
			POS = SR.split()[1]
			RSID = SR.split()[5]
			A1 = SR.split()[3]
			A2 = SR.split()[4]
			w.write(CHR + '\t' + POS + '\t' + RSID +'\t' + A1 + '\t' + A2 + '\n')
		except subprocess.CalledProcessError:
			pass
'''
SNP annotation results:
CHR	POS	RSID	A1	A2
2	17454464	rs59977443	A	T
2	105269504	rs2241797	T	A
11	17045071	rs7930058	C	T
2	82055800	rs62152638	G	A
'''
'''
# scripts for merge multi-process generated files:

fold_path = '../../datasets/middlefile/annovar_annotation/'
output_path = '../../datasets/middlefile/dbSNP_rsid_annotation.txt'

with open(output_path, 'w+') as w:
    # write file head first
    w.write('CHR\tPOS\tRSID\tA1\tA2\n')
    # merge files
    for i in range(29):
        with open(fold_path + 'dbSNP_rsid_annotation' + str(i+1) + '.txt') as r:
            line = r.readline() # ignore file head
            line = r.readline()
            while line:
                w.write(line.strip() + '\n')
                line = r.readline()

'''
#'''

# step 3 CpG and SNP annotation to methylation position of meQTLs with [CHR] splitted
#'''
meqtl_path = '../../datasets/middlefile/single_meQTL.txt'
snp_annotation_path = '../../datasets/middlefile/dbSNP_rsid_annotation.txt'
cpg_annotation_path = '../../datasets/middlefile/clean_epic/cpg_all.pkl'
save_path = '../../datasets/middlefile/meQTL_annotation_CpG_SNP_raw.pkl'

main_table = pd.read_csv(meqtl_path, sep='\t')
print(main_table.head())
'''
          CpG          SNP Ref Alt   Beta      SE         Pvalue
0  cg11913416    rs1262461   A   G -1.179  0.0348  5.600000e-119
1  cg25722041  rs114812780   G   T  1.233  0.0369  3.160000e-117
2  cg25722041   rs12568293   T   C  1.233  0.0369  3.160000e-117
3  cg25722041   rs12117910   A   C  1.233  0.0369  3.170000e-117
4  cg25722041   rs12403339   G   A  1.233  0.0369  3.170000e-117
'''
# add annotation columns
main_table['CpG_CHR'] = 'space'
main_table['CpG_POS'] = -1
main_table['SNP_CHR'] = 'space'
main_table['SNP_POS'] = -1
main_table['SNP_check_A1'] = 'space'
main_table['SNP_check_A2'] = 'space'
main_table = main_table.drop(columns=['Pvalue'])  
main_table = main_table.drop(columns=['SE']) 
print(main_table.head())

# load CpG annotation
cpg_table = pd.read_pickle(cpg_annotation_path)
indexed_cpg = cpg_table.set_index('IlmnID')
print(indexed_cpg.head())
'''
           CHR_hg38  Start_hg38   End_hg38
IlmnID
cg07881041    chr19     5236004    5236006
cg18478105    chr20    63216297   63216299
cg23229610     chr1     6781064    6781066
cg03513874     chr2   197438741  197438743
cg09835024     chrX    24054522   24054524
'''
# load SNP annotation
snp_table = pd.read_csv(snp_annotation_path,sep='\t')
indexed_snp = snp_table.set_index('RSID')
print(indexed_snp.head())
'''
            CHR        POS A1 A2
RSID
rs59977443    2   17454464  A  T
rs2241797     2  105269504  T  A
rs7930058    11   17045071  C  T
rs62152638    2   82055800  G  A
rs7298494    12   10460666  G  T
'''
for i in tqdm(range(main_table.shape[0])): # 12424038, need long time (nearly 260 hours)
    # CpG annotation
	try:
		cpg = main_table['CpG'][i]
		main_table.loc[i, 'CpG_CHR'] = indexed_cpg[indexed_cpg.index==cpg]['CHR_hg38'].values[0][3:]
		main_table.loc[i, 'CpG_POS'] = indexed_cpg[indexed_cpg.index==cpg]['Start_hg38'].values[0]
	except IndexError:
		continue        
    # SNP annotation
	snp = main_table['SNP'][i]
	try:
		main_table.loc[i, 'SNP_CHR'] = indexed_snp[indexed_snp.index==snp]['CHR'].values[0]
		main_table.loc[i, 'SNP_POS'] = indexed_snp[indexed_snp.index==snp]['POS'].values[0]
		main_table.loc[i, 'SNP_check_A1'] = indexed_snp[indexed_snp.index==snp]['A1'].values[0]
		main_table.loc[i, 'SNP_check_A2'] = indexed_snp[indexed_snp.index==snp]['A2'].values[0]
	except IndexError:
		continue
print(main_table.head())
main_table.to_pickle(save_path)
'''
         CpG          SNP Ref Alt   Beta CpG_CHR   CpG_POS SNP_CHR  SNP_POS SNP_check_A1 SNP_check_A2
4  cg25722041   rs12403339   G   A  1.233       1   8563413       1  8498232            G            A
'''

#'''

# step 4 remove wrong items
#'''


#'''