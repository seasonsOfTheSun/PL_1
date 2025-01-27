mkdir data/
mkdir data/raw/
cd data/raw/

echo $pwd
# Download the STRING PPI
wget https://stringdb-static.org/download/protein.physical.links.v11.5/9606.protein.physical.links.v11.5.txt.gz
gunzip 9606.protein.physical.links.v11.5.txt.gz

# Download info to translate ensembl to entrez ids.
# sometime in the futuure 113 might have to be incremented to 114
wget https://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.113.entrez.tsv.gz
gunzip Homo_sapiens.GRCh38.113.entrez.tsv.gz


# you can uuse this code to convert e.g. gene symbols to entrez and vice versa
wget https://ftp.ebi.ac.uk/pub/databases/genenames/out_of_date_hgnc/json/locus_groups/protein-coding_gene.json
# use the name mapper to convert whatever naming system we have
# in the gene list can be converted to entrez
