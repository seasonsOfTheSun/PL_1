

# for the time being we'll get our disease gene associations 
# from https://diseases.jensenlab.org/Downloads
# specifically the filtered knowledge channel


from ppi_tools import make_mapping_dict
import pandas as pd

def test():
    


    disease_db = load_from_jensen()
    
    diseases = disease_db.list_diseases()
    disease  = "Cancer"
    disease_db.get_genes(disease)
    


class DiseaseDataBase:
    
    def __init__(self,database_name,gene_identifier = "symbol"):
        
        self.database = database_name
        self.gene_identifier = gene_identifier
        
        self.gene_id_maps = {}
    
    
    def get_genes(self,disease):
        raw_genes = self.get_raw_genes(disease)
        try:
            gene_map = self.gene_id_maps[(self.raw_gene_identifier,self.gene_identifier)]
        except KeyError:  
            gene_map = make_mapping_dict(self.raw_gene_identifier,self.gene_identifier)
            self.gene_id_maps[(self.raw_gene_identifier,self.gene_identifier)] = gene_map
            
        out = []
        for gene in raw_genes:
            
            try:
                out.append(gene_map[gene])
            except KeyError as e:
                pass
        return out
            
            
        
        

    
class DiseaseFromDataFrame(DiseaseDataBase):
    
    def __init__(self,df,
                 database_name,
                 disease_identifier_col,
                 gene_identifier_col):
        super().__init__(database_name)
        self.df = df
        self.dic = disease_identifier_col
        self.gic = gene_identifier_col
        
        #self.raw_gene_identifier = "ensembl_protein_id"
        self.raw_gene_identifier = "symbol"

        
    def list_diseases(self):
        #return list(self.df[self.dic].unique())
        return self.df.groupby(self.dic).count()[0].sort_values(ascending=False)
        
        
    def get_raw_genes(self,disease):
        
        select = (self.df[self.dic] == disease)
        return self.df.loc[select,self.gic]
    
    
def load_from_jensen():
    df = pd.read_csv("data/raw/human_disease_knowledge_filtered.tsv",
                 sep="\t",
                 header=None)

    disease_db = DiseaseFromDataFrame(df,
                         "jensen_knowledge_filtered",
                         3,1)
    disease_db.raw_gene_identifier = "symbol"
    
    return disease_db

test()