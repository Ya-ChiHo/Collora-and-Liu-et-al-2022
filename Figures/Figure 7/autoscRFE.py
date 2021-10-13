from anndata import read_h5ad
import pandas as pd
from scRFE.scRFE import scRFE
import scanpy as sc 
from random import randint
import argparse
#goal of this script is to wrap several replicates of scRFE at once, writing each to a file, it takes in the number of replicates, the anndata file, and the output directory 
parser = argparse.ArgumentParser(description='RunsscRFE')
parser.add_argument('--iterations',"-I", metavar='I', type=int, nargs=1,
                    help='an integer for how many times to repeat process')
parser.add_argument('--datafile',"-D", metavar='D', type=str, nargs=1,
                    help='file to read in')
parser.add_argument('--outdir',"-O", metavar='O', type=str, nargs=1,
                    help='file to print out')

args = parser.parse_args()

datafile=args.datafile[0]
number=args.iterations[0]
outputfile=args.outdir[0]

adata = read_h5ad(datafile)  # read in adata

# basic filtering (optional)
sc.pp.filter_cells(adata, min_genes=250)
sc.pp.filter_genes(adata, min_cells=3)
#run the actual estimator 
i=0
while(i<number):
    rand=randint(1, 1e9)
    adata.obs['HIVclone']=[str(i) for i in list(adata.obs_vector('HIVclone'))]
    topFeaturesDF, score = scRFE(adata =  adata, classOfInterest = "HIVclone", randomState = rand)
    #save the results
    topFeaturesDF[['True','True_gini']].to_csv(outputfile+"/"+str(rand)+".csv", index=False)
    i+=1

