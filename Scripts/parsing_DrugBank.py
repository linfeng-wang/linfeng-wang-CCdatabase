"""
TITLE :parsing_DrugBank
DESCRIPTION : Parsing full database of DrugBank to extract information
                about drugs that affects a specific organism.
                Information retrieved:
                    - name
                    - synonyms
                    - classification (kingdom and superfamily)
                    - drug-interactions (with other drugs)
                    - external-identifiers (to connect to other sources)
                    - pathways
                    - targets (if polypeptides)
                        - target_name
                        - target_uniprot
                        - target_gene_name
                        - action (of the drug over the target)
                        - cell_loc (cell localitation)
            To get the tree structure of the xml file,
                see drugBank_tree_structure.txt
"""
#%%
# Classes #
class Drug:
    """
    docstring for Drug.
    """
    def __init__(self, features):

        self.id = features['id']
        self.name = features['name']
        #self.synonyms = features['synm']
        #self.kingdom = features['kgd']
        #self.superclass = features['sclass']
        self.interaction = features['itrc']
        #self.external_id = features['ext_id']
        self.pathways = features['pathways']
        self.mechanism_of_action = features['mechanism_of_action']
        self.metabolism = features['metabolism']
        self.pharmacodynamics = features['pharmacodynamics']
        self.absorption = features['absorption']
        self.protein_binding = features['protein_binding']
        self.atc_codes = features['atc_codes']
        self.target = []
        self.enzyme_list = features['enzyme_list']
        self.carrier_list = features['carrier_list']
        self.transporter_list = features['transporter_list']

    def getDrugfeatures(self):
        drug_dict = {"dg_id":self.id,
                    "dg_name":self.name,
                    #"dg_synm":self.synonyms,
                    #"dg_kingdom":self.kingdom,
                    #"dg_superclass":self.superclass,
                    "dg_interactions":self.interaction,
                    #"dg_ext_id":self.external_id,
                    "dg_pathways":self.pathways,
                    "dg_mechanism_of_action": self.mechanism_of_action,
                    "dg_metabolism" : self.metabolism,
                    "dg_pharmacodynamics" : self.pharmacodynamics,
                    "dg_absorption" : self.absorption,
                    "dg_protein_binding" : self.protein_binding,
                    "dg_atc_codes" : self.atc_codes,
                    "dg_enzyme_list" : self.enzyme_list,
                    "dg_carrier_list" : self.carrier_list,
                    "dg_transporter_list" : self.transporter_list
                    }
        
        return drug_dict

    def addTarget(self, feature_target):
        self.target.append(feature_target)

# Parameters and required variables #

dB_file = '../DB/drugbank_database.xml'
organism = 'Humans'
saveFile = '../DB/drugbank_extracted.csv'

# Main script #
'''
Get targets from drugBank database for the drugs on the ic50 file

'''
#%%
import xml.etree.ElementTree as ET
import time
from tqdm import tqdm
import pandas as pd
#%%

xtree = ET.parse(dB_file)
xroot = xtree.getroot()
drugs = list(xroot)

#%%
drug_targets = []
for i in tqdm(range(len(drugs))):
    drug = drugs[i]
    idDB = drug[0].text # Drug Bank ID

    for idx,feature in enumerate(drug):
        if 'name' in str(feature): # drug name
            drug_name = drug[idx].text

        # if 'synonyms' in str(feature): # drug's synonyms
        #     drug_synm = ';'.join([synm.text \
        #                             for synm in list(drug[idx])])

        # if 'classification' in str(feature): #type of drug
        #     drug_class_kingdom = list(drug[idx])[2].text
        #     drug_class_superclass = list(drug[idx])[3].text

        if 'drug-interactions' in str(feature): #interaction other drugs
            drug_interaction = ';'.join([di[0].text for di in list(drug[idx])])

        if 'external-identifiers' in str(feature): #other drug's IDs
            aux = [ext_id[0].text + ":" + ext_id[1].text \
                                        for ext_id in list(drug[idx])]
            drug_external_id = ';'.join(aux)
        if 'pathway' in str(feature): #related pathway
            drug_pathway = ';'.join([pathway[1].text for pathway in list(drug[idx])])
        
        
        if 'mechanism-of-action' in str(feature): 
            drug_mechanism_of_action = drug[idx].text
            #print(drug_mechanism_of_action)
            
        if 'metabolism' in str(feature): 
            drug_metabolism = drug[idx].text
            #print(metabolism)

        if 'pharmacodynamics' in str(feature): 
            drug_pharmacodynamics = drug[idx].text
            #print(drug_pharmacodynamics)
        
        if 'absorption' in str(feature): 
            drug_absorption = drug[idx].text
            #print(drug_pharmacodynamics)
        
        if 'protein-binding' in str(feature): 
            drug_protein_binding = drug[idx].text
            #print(drug_pharmacodynamics)        
            
        if 'atc-codes' in str(feature): 
            if len([di.attrib for di in list(drug[idx])]) > 0:
                codes = [di.attrib for di in list(drug[idx])]
                code_list = []
                for x in codes:
                    code_list.append(x['code'])
                drug_atc_codes = ';'.join([c for c in code_list])
            else:
                drug_atc_codes = None
            #print(drug_atc_codes)        

        if 'targets' in str(feature): #if polypeptide, drug's targets
            targets = list(drug[idx])
            
        if 'enzymes' in str(feature): 
            enzyme_list = []
            enzymes = list(drug[idx])
            for enzyme in enzymes:
                for ind, x in enumerate(enzyme):
                    if "name" in str(x):
                        enzyme_list.append(x.text)
            enzyme_list = ';'.join([c for c in enzyme_list])
                
        if 'carriers' in str(feature): 
            carriers= list(drug[idx])
            carrier_list = []
            carriers = list(drug[idx])
            for carriers in carriers:
                for ind, x in enumerate(carriers):
                    if "name" in str(x):
                        carrier_list.append(x.text)
            carrier_list = ';'.join([c for c in carrier_list])
                            
        if 'transporters' in str(feature): 
            transporters = list(drug[idx])
            transporter_list = []
            transporters = list(drug[idx])
            for transporters in transporters:
                for ind, x in enumerate(transporters):
                    if "name" in str(x):
                        transporter_list.append(x.text)
            transporter_list = ';'.join([c for c in transporter_list])

    # get all drug-related information in a dictionary
    drug_dict = {"id":idDB,
                "name":drug_name,
                #"synm":drug_synm,
                #"kgd":drug_class_kingdom,
                #"sclass":drug_class_superclass,
                "itrc":drug_interaction,
                #"ext_id":drug_external_id,
                "pathways": drug_pathway,
                "mechanism_of_action" : drug_mechanism_of_action,
                "metabolism" : drug_metabolism,
                "pharmacodynamics" : drug_pharmacodynamics,
                "absorption" : drug_absorption,
                "protein_binding" : drug_protein_binding,
                "atc_codes" : drug_atc_codes,
                "enzyme_list" : enzyme_list,
                "carriers_list": carrier_list,
                "transporters_list" : transporter_list
                }
    
    drug = Drug(drug_dict)

    # get information of polypeptide targets
    if len(targets) > 0:
        for target in targets:
            idx_pep = None
            # get indexes
            for idx,feature in enumerate(target): # check features of targets
                if 'organism' in str(feature):
                    idx_org = idx
                if 'name' in str(feature):
                    idx_name = idx
                if 'actions' in str(feature):
                    idx_act = idx
                if 'polypeptide' in str(feature):
                    idx_pep = idx

            # Get information for polypeptide
            if target[idx_org].text == organism:

                target_name = target[idx_name].text

                actions = ';'.join([action.text
                                    for action in list(target[idx_act])])

                # Get information for polypeptide
                if idx_pep is not None: #if there is polypeptide's info...
                    for idx,feature in enumerate(target[idx_pep]):
                        if 'gene-name' in str(feature):
                            gene_name = target[idx_pep][idx].text
                        if 'cellular-location' in str(feature):
                            cell_loc = target[idx_pep][idx].text
                        if 'external-identifiers' in str(feature):
                            for ext_id in list(target[idx_pep][idx]):
                                if ext_id[0].text == "UniProtKB":
                                    uniprot = ext_id[1].text
                else:
                    gene_name = None
                    action = None
                    cell_loc = None
                    uniprot = None

                row = {
                        "dg_id":drug.id,
                        "dg_name":drug.name,
                        #"dg_synm":drug.synonyms,
                        #"dg_kingdom":drug.kingdom,
                        #"dg_superclass":drug.superclass,
                        "dg_interactions":drug.interaction,
                        "dg_ext_id":drug.external_id,
                        "dg_pathways":drug.pathways,
                        "target_name":target_name,
                        "target_uniprot":uniprot,
                        "target_gene_name":gene_name,
                        "action":actions,
                        "cell_loc":cell_loc,
                        "mechanism_of_action" : drug.mechanism_of_action,
                        "metabolism" : drug.metabolism,
                        "pharmacodynamics" : drug.pharmacodynamics,
                        "absorption" : drug.absorption,
                        "protein_binding" : drug.protein_binding,
                        "atc_codes" : drug.atc_codes,
                        "enzymes" : drug.enzyme_list,
                        "carriers" : drug.carrier_list,
                        "transporters" : drug.transporter_list
                        }

                drug_targets.append(row)


dt = pd.DataFrame.from_dict(drug_targets, orient='columns')
dt.shape
dt.to_csv(saveFile)

# %%
