import sys, getopt
import pandas as pd
import numpy as np
import glob
import os 
import scanpy as sc

def main(argv):
   data_h5 = ''
   solo_outdir = ''
   output_tsv = ''
   try:
      opts, args = getopt.getopt(argv,"hd:s:o:",["data_h5=","solo_out_dir=","output_tsv="])
   except getopt.GetoptError:
      print('test.py -d <data_h5> -s <solo_out_dir> -o <output_tsv>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test.py -d <data_h5> -s <solo_out_dir> -o <output_tsv>')
         sys.exit()
      elif opt in ("-d", "--data_h5"):
         data_h5 = arg
      elif opt in ("-s", "--solo_outdir"):
         solo_outdir = arg
      elif opt in ("-o", "--ofile"):
         output_tsv = arg
   print('Input file data h5 is ', data_h5)
   print('Input file solo_outdir is ', solo_outdir)
   print('Output file is ', output_tsv)


   data_10x_h5 = sc.read_10x_h5(data_h5)
   data_10x_h5.var_names_make_unique()
   #data_10x_h5.write(work_dir + 'filtered_feature_bc_matrix2.h5ad') 
   data = {'cell_barcode':data_10x_h5.obs_names,
           'is_doublet':np.load(solo_outdir + '/' + 'is_doublet.npy'),
           'preds':np.load(solo_outdir + '/' + 'preds.npy'),
           'softmax_scores':np.load(solo_outdir + '/' + 'softmax_scores.npy'),
           'logit_scores':np.load(solo_outdir + '/' + 'logit_scores.npy')} 
  
   # Create DataFrame 
   df = pd.DataFrame(data) 
   df.to_csv(output_tsv, sep='\t', index=False)

if __name__ == "__main__":
   main(sys.argv[1:])
