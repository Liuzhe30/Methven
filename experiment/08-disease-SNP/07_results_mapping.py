# mapping results
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

snp_list = ['rs2476601','rs3806624', 'rs7731626', 'rs2234067','rs2233424','rs947474','rs3824660','rs968567','rs3218251']
pred_results_path = '../../model/pred_results_disease_SNP/'
full_t0_path = 'data/final_t0/'
full_t24_path = 'data/final_t24/'
output_path = 'data/pred_merged/'

model_size = {'small':10_000,'large':100_000}
model_size = {'small':10_000}

for snp in snp_list:
    for model in model_size.keys():
        full_t0_data = pd.read_pickle(full_t0_path + snp + '_' + model + '.dataset')
        full_t24_data = pd.read_pickle(full_t24_path + snp + '_' + model + '.dataset')
        t0_npy = np.load(pred_results_path + snp + '_' + model + '_predict_t0.npy')
        t24_npy = np.load(pred_results_path + snp + '_' + model + '_predict_t24.npy')

        full_t0_data['t0'] = 0.0
        full_t24_data['t24'] = 0.0
        for i in range(len(full_t0_data)):
            full_t0_data['t0'][i] = t0_npy[i][0]
        for i in range(len(full_t24_data)):
            full_t24_data['t24'][i] = t24_npy[i][0]

        # crossmap
        cpg_list1 = full_t0_data['CpG'].tolist()
        cpg_list2 = full_t24_data['CpG'].tolist()
        final_list = [x for x in cpg_list1 if x in cpg_list2]

        final_df = pd.DataFrame()
        for item in final_list:
            cpg = full_t0_data[full_t0_data['CpG']==item]['CpG'].values[0]
            snp = full_t0_data[full_t0_data['CpG']==item]['SNP'].values[0]
            ref = full_t0_data[full_t0_data['CpG']==item]['Ref'].values[0]
            alt = full_t0_data[full_t0_data['CpG']==item]['Alt'].values[0]
            chr = full_t0_data[full_t0_data['CpG']==item]['CHR'].values[0]
            cpg_pos = full_t0_data[full_t0_data['CpG']==item]['CpG_POS'].values[0]
            snp_pos = full_t0_data[full_t0_data['CpG']==item]['SNP_POS'].values[0]
            t0 = full_t0_data[full_t0_data['CpG']==item]['t0'].values[0]
            t24 = full_t24_data[full_t24_data['CpG']==item]['t24'].values[0]
            delta = t24 - t0
            final_df = final_df._append({'CpG':cpg,'SNP':snp,'Ref':ref,'Alt':alt,'CHR':chr,'CpG_POS':cpg_pos,
                                            'SNP_POS':snp_pos,'t0':t0,'t24':t24,'delta':delta},ignore_index=True)
        print(final_df.head())
        final_df.to_csv(output_path + snp + '_' + model + '.csv')

'''
          CpG        SNP Ref Alt CHR   CpG_POS   SNP_POS        t0       t24  \
0  cg00660737  rs3218251   T   A  22  37546043  37545505 -2.064718 -1.780150
1  cg07278285  rs3218251   T   A  22  37554775  37545505 -0.111967  0.490153
2  cg11619069  rs3218251   T   A  22  37543486  37545505  1.022392  0.692455
3  cg18060391  rs3218251   T   A  22  37536870  37545505 -3.512755 -2.801362
4  cg00154335  rs3218251   T   A  22  37554221  37545505  1.773641  0.025784

      delta
0  0.284569
1  0.602120
2 -0.329937
3  0.711393
4 -1.747857
'''