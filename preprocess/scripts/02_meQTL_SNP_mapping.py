import os
import pandas as pd
import numpy as np
import subprocess

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
			RSID = SR.split()[2]
			A1 = SR.split()[3]
			A2 = SR.split()[4]
			w.write(CHR + '\t' + POS + '\t' + RSID +'\t' + A1 + '\t' + A2 + '\n')
		except subprocess.CalledProcessError:
			pass
'''
meQTL annotation results:
CHR	POS	RSID	A1	A2
1	170193059	170193059	T	C
9	133000286	133000286	G	A
4	53477740	53477740	A	T
17	47403691	47403691	C	T
6	31790463	31790463	G	A
11	117570906	117570906	C	G
'''
'''
# scripts for merge multi-process generated files:

fold_path = '../../datasets/middlefile/annotation/'
output_path = '../../datasets/middlefile/dbSNP_rsid_annotation.txt'

with open(output_path, 'w+') as w:
    # write file head first
    w.write('CHR\tPOS\tRSID\tA1\tA2\n')
    # merge files
    for i in range(10):
        with open(fold_path + 'dbSNP_rsid_annotation' + str(i+1) + '.txt') as r:
            line = r.readline() # ignore file head
            line = r.readline()
            while line:
                w.write(line.strip() + '\n')
                line = r.readline()

'''
#'''

# step 3 SNP annotation to methylation position of meQTLs
#'''
snp_annotation_path = '../../datasets/middlefile/dbSNP_rsid_annotation.txt'


#'''

# step 4 split by [CHR] and sort by [POS]
#'''


#'''