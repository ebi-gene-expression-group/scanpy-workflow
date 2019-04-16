import json
import os
import gzip

class ID2gene_converter:
    def __init__(self, gtfFile):
        filelocation="id2gene.json"

        #Only parse the gtf file if a json file has not been created yet.
        if os.path.isfile(filelocation):
            with open(filelocation,"r") as infile:
                self.id_to_gene=json.load(infile)
        else:
            self.id_to_gene={}
            with gzip.open(gtfFile,'rt') as gtf:
                for line in gtf:
                    if line[0]!="#":
                        gene_name=""
                        gene_id=""
                        info= line.split("\t")[8].split(";")
                        for x in info:
                            if x[0:7]=="gene_id":
                                gene_id = x.split("\"")[1]
                            if x[1:10]=="gene_name":
                                gene_name = x.split("\"")[1]
                        self.id_to_gene[gene_id]=gene_name
            with open(filelocation, 'w') as outfile:  
                json.dump(self.id_to_gene, outfile)
    def id2gene(self,id):
        return self.id_to_gene[id]

