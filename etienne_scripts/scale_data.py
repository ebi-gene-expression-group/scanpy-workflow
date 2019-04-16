import numpy as np
import scanpy as sc
import sys

if len(sys.argv)<3:
    print("Not enough arguments. Usage: python scale_data.py infile outfile.")
infile = sys.argv[1]
outfile = sys.argv[2]
adata = sc.read(infile)
adata.X=np.log1p(adata.X)
adata.write(outfile)