#!/usr/bin/env python

from argparse import ArgumentParser
from os.path import isfile, isdir, join
from sys import exit
import numpy as np
#from pleiopred import prior_generating, coord_trimmed, pre_sumstats
from pleiopred import PleioPriors, coord_trimmed, pre_sumstats
from pleiopred import pleiopred_main, LD_PyWrapper
#import pred_main_bi_rho

# Create the master argparser and returns the argparser object
def get_argparser():
  parser = ArgumentParser(prog="PleioPred", 
                          description="Genetic Risk Prediction by joint modeling of multiple diseases and functional annotations.")
  ## Input Files
  #################### 
  # GWAS sumstats
  parser.add_argument('--sumstats_D1', required=True, help="GWAS summary stats of Disease 1")
  parser.add_argument('--sumstats_D2', required=True, help="GWAS summary stats of Disease 2")
  # Reference Genotype
  parser.add_argument('--ref_gt', required=True, 
                      help="Reference genotype, plink bed format")
  # Validation Genotype
  parser.add_argument('--val_gt', required=True, 
                      help="Validation genotype, plink bed format")

  ## Parameters
  # For LDSC
  parser.add_argument('--N1', required=True, type=int,
                      help="Sample size of GWAS 1, for LDSC")
#  parser.add_argument('--N_case1', required=True, type=int,
#                      help="Number of cases in GWAS training of D1, for LDSC")
#  parser.add_argument('--N_ctrl1', required=True, type=int,
#                      help="Number of ctrls in GWAS training of D1, for LDSC")

  parser.add_argument('--N2', required=True, type=int,
                      help="Sample size of GWAS 2, for LDSC")
#  parser.add_argument('--N_case2', required=True, type=int,
#                      help="Number of cases in GWAS training of D2, for LDSC")
#  parser.add_argument('--N_ctrl2', required=True, type=int,
#                      help="Number of ctrls in GWAS training of D2, for LDSC")
  parser.add_argument('--annotation_flag', required=True,
                        help="Annotation flag: Tier0, Tier1, Tier2 and Tier3")
  ## Temporary file output directory
  parser.add_argument('--temp_dir', default=".",
                      help="Directory to output all temporary files."
                           " If not specified, will use the current directory.") 
  ## Parameters
#  parser.add_argument('--N1', required=True, type=int,
#                      help="Sample size of the first disease GWAS")
#  parser.add_argument('--N2', required=True, type=int,
#                      help="Sample size of the second disease GWAS")
  parser.add_argument('--rho', type=float,
                      help="Tuning parameter in (-1,1)"
                           ", the genetic correlation between diseases")
  parser.add_argument('--alpha', type=str,
                      help="hyperparameter for the prior of PV, default 1,1,1,1")
  parser.add_argument('--init_PV', type=str,
                      help="hyperparameter for the prior of PV, default [0.25, 0.25, 0.25, 0.25]")
#  parser.add_argument('--init_betas', required=True,
#                      help="path to initial values (AnnoPred-inf scores), if not exist will be generated")
  parser.add_argument('--zero_jump_prob', type=float, default=0.0,
                      help="shrinkage level, default 0")
  parser.add_argument('--num_iter', type=int, default=250, 
                      help="Number of iterations for MCMC, default to 250.")
  parser.add_argument('--burn_in', type=int, default=100, 
                      help="burn-in for MCMC, default to 100.")
  parser.add_argument('--local_ld_prefix', required=True,
                      help="A local LD file name prefix"
                           ", will be created if not present")
#  parser.add_argument('--hfile', required=True,
#                      help="The output path for per-SNP heritability estimation, if not exist will be generated")
  parser.add_argument('--ld_radius', type=int,
                      help="If not provided, will use the number of SNPs in" 
                           " common divided by 3000")
  parser.add_argument('--coord_D1', required=True, 
                      help="Output H5 File for coord_genotypes of D1, if not exist will be generated")
  parser.add_argument('--coord_D2', required=True, 
                      help="Output H5 File for coord_genotypes of D2, if not exist will be generated")
  parser.add_argument('--out_anno', default="PleioPred_with_anno_out",
                      help="Output filename prefix for PleioPred-anno")
  parser.add_argument('--out_ld', default="PleioPred_without_anno_out",
                      help="Output filename prefix for PleioPred")
  parser.add_argument('--user_h1', type=float,
                      help="User-provided heritability estimation for D1")
  parser.add_argument('--user_h2', type=float,
                      help="User-provided heritability estimation for D2")

  return parser

# Check if all three files for PLINK exists
def check_plink_exist(prefix):
  suffices = ['bed', 'bim', 'fam']
  result = True
  for s in suffices:
    result = result and isfile(prefix + '.' + s)
  return result

def process_args(args):
  pdict = {}
  # sumstats
  if (isfile(args.sumstats_D1)):
    pdict['sumstats_D1'] = args.sumstats_D1
  else:
    exit("sumstats_D1 file does not exists!")
  if (isfile(args.sumstats_D2)):
    pdict['sumstats_D2'] = args.sumstats_D2
  else:
    exit("sumstats_D2 file does not exists!")

  # plink formats
  if (check_plink_exist(args.ref_gt)):
    pdict['ref_gt'] = args.ref_gt
  else:
    exit("Cannot find all reference genotype plink files!")
  if (check_plink_exist(args.val_gt)):
    pdict['val_gt'] = args.val_gt
  else:
    exit("Cannot find all validation genotype plink files!")
  
  if (isdir(args.temp_dir)):
    pdict['temp_dir'] = args.temp_dir
  else:
    exit("Directory for temporary files does not exist!")

  pdict['coord_D1'] = args.coord_D1
  pdict['coord_D2'] = args.coord_D2
  pdict['N1'] = args.N1
  pdict['N2'] = args.N2
  pdict['annotation_flag'] = args.annotation_flag
#  pdict['N_case1'] = args.N_case1
#  pdict['N_case2'] = args.N_case2
#  pdict['N_ctrl1'] = args.N_ctrl1
#  pdict['N_ctrl2'] = args.N_ctrl2
#  pdict['N1'] = pdict['N_case1'] + pdict['N_ctrl1']
#  pdict['N2'] = pdict['N_case2'] + pdict['N_ctrl2']


  if (args.rho is not None):
    if (args.rho>-1 and args.rho<1):
      pdict['rho'] = args.rho
    else:
      exit("Tuning parameter needs to be in (-1,1)!")
  else:
    pdict['rho'] = 0

  if (args.alpha is not None):
    pdict['alpha'] = [float(item) for item in args.alpha.split(',')]
  else:
    pdict['alpha'] = [1.0, 1.0, 1.0, 1.0]

  if (args.init_PV is not None):
    pdict['init_PV'] = [float(item) for item in args.init_PV.split(',')]
  else:
    pdict['init_PV'] = [0.25, 0.25, 0.25, 0.25]
      
  pdict['ld_radius'] = args.ld_radius
  pdict['local_ld_prefix'] = args.local_ld_prefix
#  pdict['hfile'] = args.hfile
  pdict['out_anno'] = args.out_anno
  pdict['out_ld'] = args.out_ld
  pdict['zero_jump_prob'] = args.zero_jump_prob
  pdict['num_iter'] = args.num_iter
  pdict['burn_in'] = args.burn_in
#  pdict['init_betas'] = args.init_betas
  pdict['user_h1'] = args.user_h1
  pdict['user_h2'] = args.user_h2
  return pdict

# Returns the path to the file with name in the temp directory
def tmp(pdict, name):
  return join(pdict['temp_dir'], name)

# Returns pdict used by coord_trimmed
def pdict_coord_trimmed(pdict, dis):
  d = {}
  d['N'] = pdict['N'+dis]
  d['gf'] = pdict['ref_gt']
  d['vgf'] = pdict['val_gt']
  d['ssf'] = pdict['sumstats_D'+dis]
  d['ssf_format'] = 'BASIC'
  d['out'] = pdict['coord_D'+dis]
  d['gf_format'] = 'PLINK'
  d['skip_coordination'] = False
  d['vbim'] = None
  d['gmdir'] = None
  d['indiv_list'] = None
  return d

# Returns partially filled pdict used by Pred
def pdict_pred_partial(pdict, mode):
  d = {}
  d['coord_D1'] = pdict['coord_D1']
  d['coord_D2'] = pdict['coord_D2']
  d['ld_radius'] = pdict['ld_radius']
  d['local_ld_prefix'] = pdict['local_ld_prefix']
  d['alpha'] = pdict['alpha']
  d['init_PV'] = pdict['init_PV']
  d['num_iter'] = pdict['num_iter']
  d['burn_in'] = pdict['burn_in']
  d['N1'] = pdict['N1']
  d['N2'] = pdict['N2']
  d['ld_radius'] = pdict['ld_radius']
  d['zero_jump_prob'] = pdict['zero_jump_prob']
  d['rho'] = pdict['rho']
  d['user_h1'] = pdict['user_h1']
  d['user_h2'] = pdict['user_h2']
  if mode == 'anno':
    d['hfile'] = pdict['h2_anno']
    d['init_betas'] = pdict['init_betas_anno']
    d['out'] = pdict['out_anno']
  else:
    d['hfile'] = pdict['h2_ld']
    d['init_betas'] = pdict['init_betas_ld']
    d['out'] = pdict['out_ld']
  return d


def main(pdict):
  print(pdict)
  # Filter SNPs
  print 'Filtering Summary Stats...'
  org_sumstats1 = pdict['sumstats_D1']
  sumstats_filtered1 = tmp(pdict, "sumstats_filtered_D1.txt")
  if not isfile(sumstats_filtered1):
    pre_sumstats.get_1000G_snps(pdict['sumstats_D1'], sumstats_filtered1)
  else:
    print 'Filtered sumstats_D1 found, start coordinating genotypes...'

  org_sumstats2 = pdict['sumstats_D2']
  sumstats_filtered2 = tmp(pdict, "sumstats_filtered_D2.txt")
  
  if not isfile(sumstats_filtered2):
    pre_sumstats.get_1000G_snps(pdict['sumstats_D2'], sumstats_filtered2)
  else:
    print 'Filtered sumstats_D2 found, start coordinating genotypes...'
  
  print 'Merging Summary Stats, extracting common snps ...'
  sumstats_merged1 = tmp(pdict, "sumstats_merged_D1.txt")
  sumstats_merged2 = tmp(pdict, "sumstats_merged_D2.txt")
  if not (isfile(sumstats_merged1) and isfile(sumstats_merged2)):
    pre_sumstats.merge_sumstats(sumstats_filtered1, sumstats_filtered2, sumstats_merged1, sumstats_merged2)
    pdict['sumstats_D1'] = sumstats_merged1
    pdict['sumstats_D2'] = sumstats_merged2
  else:
    print 'Merged sumstats files found, start coordinating genotypes...'

  # Generate coord_genotypes H5 file
  print 'Coordinate summary stats and validation/reference genotype data...'
  if not isfile(pdict_coord_trimmed(pdict, '1')['out']):
    coord_trimmed.main(pdict_coord_trimmed(pdict, '1'))
  else:
    print 'Coord file for D1 already exists! Continue calculating priors...'

  if not isfile(pdict_coord_trimmed(pdict, '2')['out']):
    coord_trimmed.main(pdict_coord_trimmed(pdict, '2'))
  else:
    print 'Coord file for D2 already exists! Continue calculating priors...'

  ldsc_result = tmp(pdict, pdict['annotation_flag'])
  #ldsc_result1 = tmp(pdict,'ldsc1.results')
  #ldsc_result2 = tmp(pdict,'ldsc2.results')
  if not isfile(ldsc_result+'_ldsc1.results'):
    LD_PyWrapper.callLDSC(
        org_sumstats1, pdict['N1'], ldsc_result+'_ldsc1', pdict['annotation_flag'])
  else:
    print 'LDSC results for D1 found! Continue calculating priors ...'

  if not isfile(ldsc_result+'_ldsc2.results'):
    LD_PyWrapper.callLDSC(
        org_sumstats2, pdict['N2'], ldsc_result+'_ldsc2', pdict['annotation_flag'])
  else:
    print 'LDSC results for D2 found! Continue calculating priors ...'

  pdict['h2_anno'] = tmp(pdict, pdict['annotation_flag']+"_h2_anno.txt")
  pdict['h2_ld'] = tmp(pdict, "h2_ld.txt")
  pdict['init_betas_anno'] = tmp(pdict, "init_betas_anno")
  pdict['init_betas_ld'] = tmp(pdict, "init_betas_ld")
  ld_r = PleioPriors.generate_prior_bi(
           pdict['coord_D1'], pdict['coord_D2'], 
           ldsc_result+'_ldsc1.results', ldsc_result+'_ldsc2.results',
           pdict['h2_anno'], pdict['h2_ld'], pdict['annotation_flag'])
  if pdict['ld_radius'] is None: 
    pdict['ld_radius'] = int(ld_r)
  print 'Starting PleioPred with annotations'
  pleiopred_main.main(pdict_pred_partial(pdict, 'anno'))
  print 'Starting PleioPred without annotations'
  pleiopred_main.main(pdict_pred_partial(pdict, 'ld'))



if __name__ == '__main__':
  args = get_argparser().parse_args()
  main(process_args(args))
