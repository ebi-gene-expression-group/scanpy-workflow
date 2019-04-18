import argparse
from scipy.io import mmread
import pandas as pd
import scanpy as sc
from ID2gene import ID2gene_converter

parser = argparse.ArgumentParser(description='Combines the information of several anndata files into one. Also modifies the index to use gene names instead of Ensembl IDs, and adds information from the sdrf file.')
parser.add_argument("-sdrf",type=str,help="The sdrf file of the experiment. If it is provided, it will be used to add information to the anndata matrix.")
parser.add_argument("-gtf",type=str,help="The gtf file. If it is provided, it will be used to convert Ensembl IDs to gene names.")
parser.add_argument("-matrix",type=str,help="The file where the anndata matrix is stored. It must at least contains the clusters and the umap coordinates.")
parser.add_argument("-tsne",type=str,help="A file containing an anndata matrix, which contains the t-SNE coordinates. If it is provided, the t-SNE coordinates will be added to the main matrix.")
parser.add_argument("-o",type=str,default="anndata_complete.h5ad",help="The name of the output file.")
args=parser.parse_args()

def add_sdrf_info(adata,sdrfpath):
    """Parse the sdrf file and add the characteristics to the anndata"""
    with open(sdrfpath) as sdrf:
        header=next(sdrf).split("\t")
        enaIndex=-1
        infoIndex={}
        infoData={}
        for i in range(len(header)):
            if header[i][:15]=="Characteristics":
                s=header[i]
                info = s[s.find("[")+1:s.find("]")] #select the substring between the brackets
                infoIndex[info]=i
                infoData[info]={}
            if header[i]=="Comment[ENA_RUN]" or header[i] == "Comment [ENA_RUN]":
                enaIndex=i
        for line in sdrf:
            linesplit=line.split("\t")
            ena=linesplit[enaIndex]
            if ena in adata.obs.index:
                for info in infoIndex:
                    infoData[info][ena]=linesplit[infoIndex[info]]

        #Quantitative information is converted to numeric values, so that cellxgene can use it as continuous metadata, instead of categorical.
        quantitative_info=["age","body mass index"]
        for info in infoIndex:
            if not info in quantitative_info:
                adata.obs[info]=pd.Series(infoData[info])
            else:
                adata.obs[info]=pd.to_numeric(pd.Series(infoData[info]))
            

if __name__=='__main__':
    adata=sc.read(args.matrix)
    if args.tsne is not None:
        adataTSNE = sc.read(args.tsne)
        adata.obsm["X_tsne"]=adataTSNE.obsm["X_tsne"]


    #adata.obs["louvain"] = adata.obs["louvain_r0.1"]  #cellbrowser needs an observation "louvain", but cellxgene does not.
    #adataMarkers = sc.read(folder+"counts_mtx_clusters_0.3_markers.h5ad")
    #adata.uns["rank_genes_groups"] = adataMarkers.uns["rank_genes_groups"]

    #Use gene name instead of Ensembl IDs.
    if args.gtf is not None:
        converter=ID2gene_converter(args.gtf)
        adata.raw.var.index = adata.raw.var.index.map(lambda x: converter.id2gene(x))
        adata.var.index = adata.var.index.map(lambda x: converter.id2gene(x))

    #Add some information from the sdrf file (disease, individual, age...)
    add_sdrf_info(adata,args.sdrf)

    #Create the output file
    adata.write(args.o)