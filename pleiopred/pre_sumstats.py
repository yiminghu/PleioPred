import numpy as np
import h5py

def get_1000G_snps(sumstats, out_file):
    sf = np.loadtxt(sumstats,dtype=str,skiprows=1)
    h5f = h5py.File('ref/Misc/1000G_SNP_info.h5','r')
    rf = h5f['snp_chr'][:]
    h5f.close()
    ind1 = np.in1d(sf[:,1],rf[:,2])
    ind2 = np.in1d(rf[:,2],sf[:,1])
    sf1 = sf[ind1]
    rf1 = rf[ind2]
    ### check order ###
    if sum(sf1[:,1]==rf1[:,2])==len(rf1[:,2]):
        print 'Good!'
    else:
        print 'Shit happens, sorting sf1 to have the same order as rf1'
        O1 = np.argsort(sf1[:,1])
        O2 = np.argsort(rf1[:,2])
        O3 = np.argsort(O2)
        sf1 = sf1[O1][O3]
    out = ['hg19chrc snpid a1 a2 bp or p'+'\n']
    for i in range(len(sf1[:,1])):
        out.append(sf1[:,0][i]+' '+sf1[:,1][i]+' '+sf1[:,2][i]+' '+sf1[:,3][i]+' '+rf1[:,1][i]+' '+sf1[:,5][i]+' '+sf1[:,6][i]+'\n')
    ff = open(out_file,"w")
    ff.writelines(out)
    ff.close()


def merge_sumstats(ss1, ss2, output1, output2):
	print 'Merging training summary stats of two diseases, extracting overlapping snps!'
	sumstats1 = np.loadtxt(ss1,dtype=str,skiprows=1)
	sumstats2 = np.loadtxt(ss2,dtype=str,skiprows=1)
	ovp = np.intersect1d(sumstats1[:,1], sumstats2[:,1])
	st1 = sumstats1[np.in1d(sumstats1[:,1], ovp)]
	st2 = sumstats2[np.in1d(sumstats2[:,1], ovp)]
	if np.sum(st1[:,1]==st2[:,1])!=len(ovp):
		print 'Shit happens, sorting ss1 to have the same order as ss2'
		O1 = np.argsort(st1[:,1])
		O2 = np.argsort(st2[:,2])
		O3 = np.argsort(O2)
		st1 = st1[O1][O3]
	else:
		print 'Good!'

	shit = np.logical_not(np.logical_or(np.logical_and(st1[:,2]==st2[:,2],st1[:,3]==st2[:,3]),np.logical_and(st1[:,2]==st2[:,3],st1[:,3]==st2[:,2])))
	print 'Number of snps with flipped alleles: ' + str(np.sum(shit))
	prb2_1 = st1[shit,2]
	prb2_2 = st1[shit,3]

	flipped1 = np.empty(len(prb2_1), dtype=str)
	flipped2 = np.empty(len(prb2_2), dtype=str)
	for i in range(len(prb2_1)):
		flipped1[i] = flip(prb2_1[i])
		flipped2[i] = flip(prb2_2[i])

	st1[shit,2] = flipped1
	st1[shit,3] = flipped2
	
	same = np.logical_and(st1[:,2]==st2[:,2], st1[:,3]==st2[:,3])
	reverse = np.logical_and(st1[:,2]==st2[:,3], st1[:,3]==st2[:,2])
	st1_2 = st1[np.logical_or(same, reverse)]
	st2_2 = st2[np.logical_or(same, reverse)]
	reverse1 = np.logical_and(st2_2[:,2]==st1_2[:,3], st2_2[:,3]==st1_2[:,2])
	tmp_a1 = st2_2[reverse1,2]
	tmp_a2 = st2_2[reverse1,3]
	tmp_or = st2_2[reverse1,5]
	st2_2[reverse1,2] = tmp_a2
	st2_2[reverse1,3] = tmp_a1
	st2_2[reverse1,5] = 1.0/tmp_or.astype(float)
	
	print 'Write merged files to '+output1+' and '+output2
	out1 = ['hg19chrc snpid a1 a2 bp or p'+'\n']
	for i in range(len(st1_2[:,1])):
		out1.append(st1_2[:,0][i]+' '+st1_2[:,1][i]+' '+st1_2[:,2][i]+' '+st1_2[:,3][i]+' '+st1_2[:,4][i]+' '+st1_2[:,5][i]+' '+st1_2[:,6][i]+'\n')
	ff = open(output1,"w")
	ff.writelines(out1)
	ff.close()
	out2 = ['hg19chrc snpid a1 a2 bp or p'+'\n']
	for i in range(len(st2_2[:,1])):
		out2.append(st2_2[:,0][i]+' '+st2_2[:,1][i]+' '+st2_2[:,2][i]+' '+st2_2[:,3][i]+' '+st2_2[:,4][i]+' '+st2_2[:,5][i]+' '+st2_2[:,6][i]+'\n')
	ff = open(output2,"w")
	ff.writelines(out2)
	ff.close()

def flip(allele):
	dc = {'A':'T','T':'A','C':'G','G':'C'}
	if len(allele)==1:
		return dc[allele]
	else:
		return 'Err'

