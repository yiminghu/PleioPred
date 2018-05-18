import pandas as pd
import numpy as np
import argparse
'''
1. exclude ambiguous SNP
2. get the overlap of three SNP sets
3. determine which SNP needs to flip sign
'''

ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
valid_nts = set(['A','T','C','G'])
def merge_ss_bim(ss1_rs_allele, ss2_rs_allele, bim_file, temp_path):
	ss1 = pd.read_table(ss1_rs_allele, sep=' ')
	ss1.rename(columns = {'snpid':'SNP','a1':'A1_ss1','a2':'A2_ss1'}, inplace = True)
	ss1['index'] = range(ss1.shape[0])
	ss2 = pd.read_table(ss2_rs_allele)
	ss2.rename(columns = {'snpid':'SNP','a1':'A1_ss2','a2':'A2_ss2'}, inplace = True)
	ss2['index'] = range(ss2.shape[0])
	bim = pd.read_table(bim_file, names = ["CHROM", "SNP", "GP", "POS", "A1", "A2"])
	bim['index'] = range(bim.shape[0])
	## merge to get the sorted overlapped snp list
	tmp = pd.merge(ss1, ss2, how='inner', on=['SNP'])
	tmp = pd.merge(tmp, bim, how='inner', on=['SNP'])
	## generate two lists: (rs, chr, pos, ind1, ind2, ind3), 1: keep and change of sign, 0: keep and no change of sign, 
	## -1: delete, could be either in ambiguous or shit happens when at least one of ss1 and ss2 does not match bim
	important_list, get_index1, get_index2 = [], [], []
	for i in range(tmp.shape[0]):
		g_nt = [tmp['A1'].iloc[i],tmp['A2'].iloc[i]]
		ss1_nt = [tmp['A1_ss1'].iloc[i],tmp['A2_ss1'].iloc[i]]
		ss2_nt = [tmp['A1_ss2'].iloc[i],tmp['A2_ss2'].iloc[i]]
		if tuple(g_nt) in ambig_nts:
			#important_list.append([tmp['SNP'].iloc[i], tmp['CHROM'].iloc[i], tmp['POS'].iloc[i], tmp['A1'].iloc[i], tmp['A2'].iloc[i], -1, -1, -1])
			continue
		if (not g_nt[0] in valid_nts) or (not g_nt[1] in valid_nts):
			#important_list.append([tmp['SNP'].iloc[i], tmp['CHROM'].iloc[i], tmp['POS'].iloc[i], tmp['A1'].iloc[i], tmp['A2'].iloc[i], -1, -1, -1])
			continue
		tmp_ind1, tmp_ind2 = 1, 1
		flip_nts1, flip_nts2 = False, False
		os_g_nt = np.array([opp_strand_dict[g_nt[0]], opp_strand_dict[g_nt[1]]])
		if (np.all(g_nt == ss1_nt) or np.all(os_g_nt == ss1_nt)):
			tmp_ind1 = 1
		else:
			# Opposite strand nucleotides
			flip_nts1 = (g_nt[1] == ss1_nt[0] and g_nt[0] == ss1_nt[1]) or (os_g_nt[1] == ss1_nt[0] and os_g_nt[0] == ss1_nt[1])
			if flip_nts1:
				tmp_ind1 = -1
			else:
				tmp_ind1 = 0

		if (np.all(g_nt == ss2_nt) or np.all(os_g_nt == ss2_nt)):
			tmp_ind2 = 1
		else:
			# Opposite strand nucleotides
			flip_nts2 = (g_nt[1] == ss2_nt[0] and g_nt[0] == ss2_nt[1]) or (os_g_nt[1] == ss2_nt[0] and os_g_nt[0] == ss2_nt[1])
			if flip_nts2:
				tmp_ind2 = -1
			else:
				tmp_ind2 = 0

		if (tmp_ind1 == 0) or (tmp_ind2 == 0):
			#important_list.append([tmp['SNP'].iloc[i], tmp['CHROM'].iloc[i], tmp['POS'].iloc[i], tmp['A1'].iloc[i], tmp['A2'].iloc[i], -1, -1, -1])
			continue
		else:
			important_list.append([tmp['SNP'].iloc[i], tmp['CHROM'].iloc[i], tmp['POS'].iloc[i], \
				tmp['A1'].iloc[i], tmp['A2'].iloc[i], tmp_ind1, tmp_ind2])
			get_index1.append([tmp['SNP'].iloc[i],tmp['A1_ss1'].iloc[i],tmp['A2_ss1'].iloc[i]])
			get_index2.append([tmp['SNP'].iloc[i],tmp['A1_ss2'].iloc[i],tmp['A2_ss2'].iloc[i]])
	print '%d SNPs found in overlap of summary statistics and validation data! Will be used in prediction models!'%len(important_list)
	important_df = pd.DataFrame(important_list, columns = ["SNP", "CHROM", "POS", "A1", "A2", "ind_ss1", "ind_ss2"])
	tmp = pd.merge(important_df, bim, how='left', on=["SNP", "A1", "A2"])
	get_index_df1 = pd.DataFrame(get_index1, columns = ["SNP", "A1_ss1", "A2_ss1"])
	get_index_df2 = pd.DataFrame(get_index2, columns = ["SNP", "A1_ss2", "A2_ss2"])
	get_index_df1 = pd.merge(get_index_df1, ss1, how='left', on=["SNP", "A1_ss1", "A2_ss1"])
	get_index_df2 = pd.merge(get_index_df2, ss2, how='left', on=["SNP", "A1_ss2", "A2_ss2"])
	important_df['cord_ss1'] = get_index_df1['index']
	important_df['cord_ss2'] = get_index_df2['index']
	important_df['cord_bim'] = tmp['index']
	important_df.sort_values(by=['CHROM', 'POS'], inplace=True)
	print 'SNP list written to %s!'%temp_path
	important_df.to_csv(temp_path, header=True, index=False, sep='\t')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Find overlap of SNP set in ss1, ss2 and bim.')
	parser.add_argument('--ss1_path', required=True)
	parser.add_argument('--ukb_ref', required=True)
	parser.add_argument('--val_bim', required=True)
	parser.add_argument('--output_path', required=True)
	args = parser.parse_args()
	merge_ss_bim(args.ss1_path, args.ukb_ref, args.val_bim, args.output_path)

'''
test case
ss1_rs_allele = '/gpfs/loomis/project/fas/zhao/yh367/EsmlPred/Data/train/ATH/GABRIEL_lite_N_26475.txt'
ss2_rs_allele = '/gpfs/loomis/project/fas/zhao/yh367/ukb_2419_traits_rs_allele'
bim_file = '/gpfs/loomis/project/fas/zhao/yh367/EsmlPred/Data/test/ATH/phs000788/eur_ath_gr_removed_mapped.bim'
temp_path = '/gpfs/loomis/project/fas/zhao/yh367/EsmlPred/temp/keep_flip.list'
merge_ss_bim(ss1_rs_allele, ss2_rs_allele, bim_file, temp_path)
python merge_input_ss_ukb.py \
 --ss1_path /gpfs/loomis/project/fas/zhao/yh367/EsmlPred/Data/train/ATH/GABRIEL_lite_N_26475.txt \
 --ukb_ref /gpfs/loomis/project/fas/zhao/yh367/ukb_2419_traits_rs_allele \
 --val_bim /gpfs/loomis/project/fas/zhao/yh367/EsmlPred/Data/test/ATH/phs000788/eur_ath_gr_removed_mapped.bim \
 --output_path /gpfs/loomis/project/fas/zhao/yh367/EsmlPred/temp/keep_flip.list
'''

