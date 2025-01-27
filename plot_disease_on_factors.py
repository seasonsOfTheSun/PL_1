import os
import pandas as pd
from diseases import load_from_jensen
import scipy.stats
import tqdm
import numpy as np

import matplotlib.pyplot as plt

def main(factor_dataframe_folder):
    
    disease_db = load_from_jensen()
    
    diseases = disease_db.list_diseases()

    
    factor_df = "rwd_factors"
    dataframe_filename = os.path.join(factor_dataframe_folder, f"{factor_df}.csv")
    factors_df = pd.read_csv(dataframe_filename, index_col=0)
    

    
    out = {}
    for disease in tqdm.tqdm(diseases.index):
        out_disease = {}
        for factor_column in factors_df.columns[:4]:
            
            
            genes = list(disease_db.get_genes(disease))
            #print(disease)
            #print(genes)
            
            factors = factors_df[factor_column]
            
            genes_in_network = list(set(genes) & set(factors_df.index))
            #print(genes)
            
            if len(genes_in_network) == 0:
                continue
            #fig,ax = plt.subplots()
            #ax.hist(factors.loc[genes_in_network],density=True,bins=100,log=True)
            #ax.hist(factors,density=True,bins=100,log=True)
            #fig.show()
            
            mean_ = np.mean(factors.values)
            std_ = np.std(factors.values)
            out_disease[factor_column] = (np.mean(factors[genes_in_network].values) - mean_)/std_
        out[disease] = out_disease

    df = pd.DataFrame(out)
    
    output_folder = "data/intermediate/disease_factor_z_score"
    os.makedirs(output_folder,exist_ok=True)
    filename = f"{factor_df}_jensen_disease_db"
    full_filename = os.path.join(output_folder,filename)
    df.to_csv(full_filename)
    
   
            
    
    
    
    
    
    
     

main("data/intermediate/factor_dataframes")    

