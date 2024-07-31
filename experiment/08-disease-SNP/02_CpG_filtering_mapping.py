import pandas as pd
import numpy as np

cpg_path = '../../datasets/middlefile/clean_epic/'
snp_path = 'data/filtered_SNP_hg19.csv'

snp_table = pd.read_csv(snp_path)
#print(snp_table.head())
'''
        RSID  CHR        POS
0  rs2476601    1  114377568
1  rs3806624    3   27764623
2  rs7731626    5   55444683
3  rs2234067    6   36355654
4  rs2233424    6   44233921
'''

# generate singel SNP datasets
for i in range(len(snp_table)):
    snp_small = pd.DataFrame()
    snp_large = pd.DataFrame()
    rsid = snp_table['RSID'][i]
    chr = str(snp_table['CHR'][i])
    pos = int(snp_table['POS'][i])
    ref = snp_table['Ref'][i]
    alt = snp_table['Alt'][i]
    cpg_data = pd.read_pickle(cpg_path + 'cpg_chr' + chr + '.pkl')
    for j in range(len(cpg_data)):
        cpg = cpg_data['IlmnID'][j]
        cpg_pos = cpg_data['CpG_pos'][j]
        distance = np.abs(int(cpg_pos) - int(pos))
        if(pos-10_000 < cpg_pos < pos+10_000):
            snp_small = snp_small._append({'CpG':cpg,'SNP':rsid,'Beta':0,'Ref':ref,'Alt':alt,'CHR':chr,'CpG_POS':cpg_pos,
                                            'SNP_POS':pos,'label':0,'distance':distance},ignore_index=True)
        elif(pos-100_000 < cpg_pos < pos-10_000 or pos+10_000 < cpg_pos < pos+100_000):
            snp_large = snp_large._append({'CpG':cpg,'SNP':rsid,'Beta':0,'Ref':ref,'Alt':alt,'CHR':chr,'CpG_POS':cpg_pos,
                                            'SNP_POS':pos,'label':0,'distance':distance},ignore_index=True)
    print(snp_small.head())
    print(len(snp_small))
    print(snp_large.head())
    print(len(snp_large))
    print()
    snp_small.to_pickle('data/snp_cpg/' + rsid + '_small.pkl')
    snp_large.to_pickle('data/snp_cpg/' + rsid + '_large.pkl')

'''
          CpG        SNP  Beta Ref Alt CHR    CpG_POS    SNP_POS  label  distance
0  cg12900790  rs2476601     0   A   G   1  114375675  114377568      0      1893
1  cg01439475  rs2476601     0   A   G   1  114372756  114377568      0      4812
2
          CpG        SNP  Beta Ref Alt CHR    CpG_POS    SNP_POS  label  distance
0  cg08266106  rs2476601     0   A   G   1  114430903  114377568      0     53335
1  cg05825720  rs2476601     0   A   G   1  114301557  114377568      0     76011
2  cg02066222  rs2476601     0   A   G   1  114354688  114377568      0     22880
3  cg13572289  rs2476601     0   A   G   1  114447746  114377568      0     70178
4  cg18973817  rs2476601     0   A   G   1  114471878  114377568      0     94310
135

          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg26322816  rs3806624     0   A   G   3  27760732  27764623      0      3891
1  cg26106948  rs3806624     0   A   G   3  27766284  27764623      0      1661
2  cg07613752  rs3806624     0   A   G   3  27756421  27764623      0      8202
3  cg15308664  rs3806624     0   A   G   3  27763797  27764623      0       826
4  cg03935055  rs3806624     0   A   G   3  27771905  27764623      0      7282
62
          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg14340541  rs3806624     0   A   G   3  27691539  27764623      0     73084
1  cg03363248  rs3806624     0   A   G   3  27754391  27764623      0     10232
2  cg21618784  rs3806624     0   A   G   3  27863205  27764623      0     98582
3  cg13707894  rs3806624     0   A   G   3  27674461  27764623      0     90162
4  cg00210003  rs3806624     0   A   G   3  27698966  27764623      0     65657
44

          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg05899688  rs7731626     0   G   A   5  55448527  55444683      0      3844
1  cg00391767  rs7731626     0   G   A   5  55442689  55444683      0      1994
2  cg22447380  rs7731626     0   G   A   5  55451073  55444683      0      6390
3  cg12849108  rs7731626     0   G   A   5  55452517  55444683      0      7834
4  cg21124310  rs7731626     0   G   A   5  55444106  55444683      0       577
12
          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg19678610  rs7731626     0   G   A   5  55529375  55444683      0     84692
1  cg11264935  rs7731626     0   G   A   5  55529124  55444683      0     84441
2  cg12725358  rs7731626     0   G   A   5  55529765  55444683      0     85082
3  cg09577966  rs7731626     0   G   A   5  55411971  55444683      0     32712
4  cg06750410  rs7731626     0   G   A   5  55412925  55444683      0     31758
65

          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg21220670  rs2234067     0   A   C   6  36356244  36355654      0       590
1  cg18455785  rs2234067     0   A   C   6  36354489  36355654      0      1165
2  cg13133961  rs2234067     0   A   C   6  36355488  36355654      0       166
3  cg22089736  rs2234067     0   A   C   6  36359367  36355654      0      3713
4  cg05298677  rs2234067     0   A   C   6  36355542  36355654      0       112
26
          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg11326036  rs2234067     0   A   C   6  36341704  36355654      0     13950
1  cg02959688  rs2234067     0   A   C   6  36283929  36355654      0     71725
2  cg04795450  rs2234067     0   A   C   6  36288118  36355654      0     67536
3  cg23008426  rs2234067     0   A   C   6  36308947  36355654      0     46707
4  cg24219671  rs2234067     0   A   C   6  36275520  36355654      0     80134
87

          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg08855022  rs2233424     0   C   A   6  44235210  44233921      0      1289
1  cg24652615  rs2233424     0   C   A   6  44243304  44233921      0      9383
2  cg08290850  rs2233424     0   C   A   6  44238228  44233921      0      4307
3  cg12996851  rs2233424     0   C   A   6  44238233  44233921      0      4312
4  cg17044843  rs2233424     0   C   A   6  44240377  44233921      0      6456
42
          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg16242498  rs2233424     0   C   A   6  44213193  44233921      0     20728
1  cg22674613  rs2233424     0   C   A   6  44310331  44233921      0     76410
2  cg07461996  rs2233424     0   C   A   6  44195227  44233921      0     38694
3  cg07826642  rs2233424     0   C   A   6  44201706  44233921      0     32215
4  cg05210969  rs2233424     0   C   A   6  44145950  44233921      0     87971
184

          CpG       SNP  Beta Ref Alt CHR  CpG_POS  SNP_POS  label  distance
0  cg15084823  rs947474     0   G   A  10  6392160  6390450      0      1710
1  cg18100746  rs947474     0   G   A  10  6390958  6390450      0       508
2  cg24277140  rs947474     0   G   A  10  6392062  6390450      0      1612
3  cg10718100  rs947474     0   G   A  10  6392069  6390450      0      1619
4  cg22423996  rs947474     0   G   A  10  6392132  6390450      0      1682
7
          CpG       SNP  Beta Ref Alt CHR  CpG_POS  SNP_POS  label  distance
0  cg09043941  rs947474     0   G   A  10  6297748  6390450      0     92702
1  cg10135213  rs947474     0   G   A  10  6333890  6390450      0     56560
2  cg27278668  rs947474     0   G   A  10  6452060  6390450      0     61610
3  cg16980393  rs947474     0   G   A  10  6438364  6390450      0     47914
4  cg20005705  rs947474     0   G   A  10  6297801  6390450      0     92649
56

          CpG        SNP  Beta Ref Alt CHR  CpG_POS  SNP_POS  label  distance
0  cg17489908  rs3824660     0   C   G  10  8101566  8104722      0      3156
1  cg01166071  rs3824660     0   C   G  10  8095687  8104722      0      9035
2  cg15187550  rs3824660     0   C   G  10  8096370  8104722      0      8352
3  cg11018337  rs3824660     0   C   G  10  8095495  8104722      0      9227
4  cg17566118  rs3824660     0   C   G  10  8095797  8104722      0      8925
67
          CpG        SNP  Beta Ref Alt CHR  CpG_POS  SNP_POS  label  distance
0  cg20673652  rs3824660     0   C   G  10  8126649  8104722      0     21927
1  cg02558407  rs3824660     0   C   G  10  8145887  8104722      0     41165
2  cg25988034  rs3824660     0   C   G  10  8093239  8104722      0     11483
3  cg08482531  rs3824660     0   C   G  10  8078127  8104722      0     26595
4  cg27370808  rs3824660     0   C   G  10  8079737  8104722      0     24985
96

          CpG       SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg07709195  rs968567     0   C   T  11  61586015  61595564      0      9549
1  cg08281583  rs968567     0   C   T  11  61595223  61595564      0       341
2  cg02563962  rs968567     0   C   T  11  61595550  61595564      0        14
3  cg20896974  rs968567     0   C   T  11  61595983  61595564      0       419
4  cg15454066  rs968567     0   C   T  11  61595377  61595564      0       187
32
          CpG       SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg07463765  rs968567     0   C   T  11  61631264  61595564      0     35700
1  cg09955919  rs968567     0   C   T  11  61534763  61595564      0     60801
2  cg11867718  rs968567     0   C   T  11  61647697  61595564      0     52133
3  cg25628257  rs968567     0   C   T  11  61559992  61595564      0     35572
4  cg06152049  rs968567     0   C   T  11  61640586  61595564      0     45022
231

          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg00660737  rs3218251     0   T   A  22  37546043  37545505      0       538
1  cg07278285  rs3218251     0   T   A  22  37554775  37545505      0      9270
2  cg11619069  rs3218251     0   T   A  22  37543486  37545505      0      2019
3  cg18060391  rs3218251     0   T   A  22  37536870  37545505      0      8635
4  cg00154335  rs3218251     0   T   A  22  37554221  37545505      0      8716
20
          CpG        SNP  Beta Ref Alt CHR   CpG_POS   SNP_POS  label  distance
0  cg17088007  rs3218251     0   T   A  22  37458698  37545505      0     86807
1  cg23671006  rs3218251     0   T   A  22  37578450  37545505      0     32945
2  cg15271026  rs3218251     0   T   A  22  37494173  37545505      0     51332
3  cg04631160  rs3218251     0   T   A  22  37561497  37545505      0     15992
4  cg23423423  rs3218251     0   T   A  22  37603155  37545505      0     57650
176
'''