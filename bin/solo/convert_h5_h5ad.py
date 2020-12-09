import scanpy as sc
import sys, getopt

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print('test.py -i <inputfile.h5> -o <outputfile.h5ad>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   print('Input file is "', inputfile)
   print('Output file is "', outputfile)


   adata = sc.read_10x_h5(
       inputfile
       # var_names='gene_ids',
       # make_unique=False,
       # var_names='gene_symbols'
   )
   adata.var_names_make_unique()
   adata.write(outputfile)  # , compression='gzip')

if __name__ == "__main__":
   main(sys.argv[1:])

