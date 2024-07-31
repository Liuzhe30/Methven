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
            final_df = final_df._append({'CpG':cpg,'SNP':snp,'Ref':ref,'Alt':alt,'CHR':chr,'CpG_POS':cpg_pos,
                                            'SNP_POS':snp_pos,'t0':t0,'t24':t24},ignore_index=True)
        print(final_df.head())
        final_df.to_csv(output_path + snp + '_' + model + '.csv')

'''
          CpG       SNP Ref Alt CHR  CpG_POS  SNP_POS        t0       t24
0  cg15084823  rs947474   G   A  10  6392160  6390450  5.808340  6.500612
1  cg18100746  rs947474   G   A  10  6390958  6390450  0.073516 -0.612236
2  cg24277140  rs947474   G   A  10  6392062  6390450  4.486728  4.792028
3  cg10718100  rs947474   G   A  10  6392069  6390450  4.904702  5.263696
4  cg22423996  rs947474   G   A  10  6392132  6390450  5.301204  6.127559
          CpG        SNP Ref Alt CHR  CpG_POS  SNP_POS        t0       t24
0  cg17489908  rs3824660   C   G  10  8101566  8104722 -0.893736 -0.806910
1  cg01166071  rs3824660   C   G  10  8095687  8104722  2.413714  0.571770
2  cg15187550  rs3824660   C   G  10  8096370  8104722 -3.281850 -2.574403
3  cg11018337  rs3824660   C   G  10  8095495  8104722  1.451097  0.640844
4  cg17566118  rs3824660   C   G  10  8095797  8104722  1.147711 -0.970443
'''