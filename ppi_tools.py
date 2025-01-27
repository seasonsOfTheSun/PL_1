import networkx as nx
import pandas as pd
import json

def build_PPI():
    ppi = nx.Graph()
    temp = pd.read_csv("data/raw/9606.protein.physical.links.v11.5.txt",sep=" ").values
    
    ensembl_entrez = make_ensembl_entrez()
    d = make_entrez_dictionaries()
    
    #print(try_ensembl_entrez)
    for i,j,w in temp:
        if i == j:
            continue # avoid self-loops
        
        try:
            entrez_i = ensembl_entrez[i.split(".")[1]]
            entrez_j = ensembl_entrez[j.split(".")[1]]
            
            symbol_i = d[entrez_i]["symbol"]
            symbol_j = d[entrez_j]["symbol"]
            
            ppi.add_edge(symbol_i,symbol_j,
                         weight=1.0
                         )
            
    
        except KeyError:
            pass
    
    d = make_symbol_dictionaries()
    for node in ppi.nodes():
        for i,v in d[node].items():
                 ppi.nodes()[node][i] = v

    return ppi

def max_connected_component(G):
    # subselect to only include nodes in the gcc
    gcc = max(nx.connected_components(G),key=len)
    return G.subgraph(gcc)

def make_entrez_dictionaries():
    d = json.load(open("data/raw/protein-coding_gene.json"))
    l = d["response"]["docs"]
    
    out = {}
    for entry in l:
        try:
            out[entry["entrez_id"]] = entry
        except KeyError:
            pass
    return out


def make_symbol_dictionaries():
    d = json.load(open("data/raw/protein-coding_gene.json"))
    l = d["response"]["docs"]
    
    out = {}
    for entry in l:
        try:
            out[entry["symbol"]] = entry
        except KeyError:
            pass
    return out

def make_ensembl_entrez():
    # dictionary that maps ensembl protein ids to entrez gene ids
    ensembl_entrez = pd.read_csv("data/raw/Homo_sapiens.GRCh38.113.entrez.tsv", sep="\t").set_index("protein_stable_id")["xref"].to_dict()
    del ensembl_entrez['-']
    return ensembl_entrez

def make_entrez_ensembl():
    d = make_ensembl_entrez()
    return {v:i for i,v in d.items()}
    

def make_mapping_dict(key_entry,value_entry):
    
    d = json.load(open("data/raw/protein-coding_gene.json"))
    l = d["response"]["docs"]

    ee_d = make_entrez_ensembl()
        
    out = {}
    for entry in l:
        
        try:
            if (key_entry == "ensembl_protein_id") or (value_entry == "ensembl_protein_id"):
                entry["ensembl_protein_id"] = ee_d[entry["entrez_id"]]
                out[entry[key_entry]] = entry[value_entry]

            else:
                out[entry[key_entry]] = entry[value_entry]
        except KeyError:
            pass
    return out

build_PPI()