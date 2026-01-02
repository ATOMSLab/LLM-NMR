import pubchempy as pcp
import time
import requests              
from urllib.parse import quote


def generate_smiles_string(name):
    
    if not name or name == "FAILED EXTRACTION":
        return "FAILED"
    
    #Try direct lookup with given name
    try:
        time.sleep(0.25) 
        results = pcp.get_properties('CanonicalSMILES', name, 'name')
        if results and len(results) > 0:
            smiles = results[0].get('CanonicalSMILES') or results[0].get('ConnectivitySMILES')
            if smiles:
                return smiles
    except Exception as e:
        print(f"  PubChem direct lookup failed for '{name}': {e}")
    
    #Get compound info from PubChem
    try:
        time.sleep(0.25)
        compounds = pcp.get_compounds(name, 'name')
        if compounds and len(compounds) > 0:
            compound = compounds[0]
            
            #Get SMILES directly from compound object
            if compound.canonical_smiles:
                smiles = compound.canonical_smiles
                return smiles
            
            #Try IUPAC name if available
            iupac_name = compound.iupac_name
            if iupac_name and iupac_name != name:                
                time.sleep(0.25) 
                results = pcp.get_properties('CanonicalSMILES', iupac_name, 'name')
                if results and len(results) > 0:
                    smiles = results[0].get('CanonicalSMILES')
                    if smiles:
                        return smiles
                
    except Exception as e:
        print(f"  PubChem compound lookup failed for '{name}': {e}")
    
    #Try OPSIN
    try:
        time.sleep(0.25)
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{quote(name)}.smi"
        response = requests.get(opsin_url, timeout=10)
        
        if response.status_code == 200:
            smiles = response.text.strip()
            # Validate - should be actual SMILES
            if smiles and not smiles.startswith('Could not') and '<' not in smiles:
                return smiles
    except Exception as e:
        print(f"  OPSIN lookup failed for '{name}': {e}")
    
    return "FAILED"
