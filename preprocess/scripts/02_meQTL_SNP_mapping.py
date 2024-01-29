import os
import pandas as pd
import numpy as np
import subprocess
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

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

# step 3 CpG and SNP annotation to methylation position of meQTLs
#'''
meqtl_path = '../../datasets/middlefile/single_meQTL.txt'
snp_annotation_path = '../../datasets/middlefile/dbSNP_rsid_annotation.txt'
cpg_annotation_path = '../../datasets/middlefile/clean_epic/cpg_all.pkl'
save_path = '../../datasets/middlefile/meQTL_annotation_CpG_SNP.pkl'

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
new_table = pd.DataFrame(columns=['CpG', 'SNP', 'Beta', 'Ref', 'Alt', 'CHR', 'CpG_POS', 'SNP_POS'])
for i in tqdm(range(main_table.shape[0])): # 12424038, need long time in single thread (500 items/min)
	CpG_CHR = 'space'
	SNP_CHR = 'space'
	try:
		# CpG annotation
		cpg = main_table['CpG'][i]
		CpG_CHR = indexed_cpg[indexed_cpg.index==cpg]['CHR_hg38'].values[0][3:]
		CpG_POS = indexed_cpg[indexed_cpg.index==cpg]['Start_hg38'].values[0]
		# SNP annotation
		snp = main_table['SNP'][i]
		SNP_CHR = indexed_snp[indexed_snp.index==snp]['CHR'].values[0]
		SNP_POS = indexed_snp[indexed_snp.index==snp]['POS'].values[0]
		SNP_check_A1 = indexed_snp[indexed_snp.index==snp]['A1'].values[0]
		SNP_check_A2 = indexed_snp[indexed_snp.index==snp]['A2'].values[0]
		# check items: null CpG | null SNP | wrong A1 | wrong A2 | different CHR
		if(CpG_CHR != 'space' and SNP_CHR != 'space' and main_table['Ref'][i] == SNP_check_A1 and main_table['Alt'][i] == SNP_check_A2 and str(CpG_CHR) == str(SNP_CHR)):
			new_table = new_table.append([{'CpG':cpg, 'SNP':snp, 'Beta':main_table['Beta'][i], 'Ref':main_table['Ref'][i], 
                                    'Alt':main_table['Alt'][i], 'CHR':CpG_CHR, 'CpG_POS':CpG_POS, 'SNP_POS':SNP_POS}], ignore_index=True)
	except IndexError:
		continue        
	
print(new_table.head())
print(new_table.shape)
new_table.to_pickle(save_path)
'''
          CpG         SNP   Beta Ref Alt CHR  CpG_POS  SNP_POS
0  cg25722041  rs12403339  1.233   G   A   1  8563413  8498232
1  cg25722041  rs28422051  1.233   T   C   1  8563413  8558278
2  cg25722041  rs11584261  1.233   C   A   1  8563413  8500306
3  cg25722041  rs12405049  1.233   C   T   1  8563413  8566058
4  cg25722041   rs2144463  1.233   C   A   1  8563413  8523234
'''

#'''