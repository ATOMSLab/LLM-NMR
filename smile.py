import cirpy
import time

def generate_smiles_string(name):
    try:
        # Pause for 15 seconds to avoid overwhelming the resolver server 
        # (recommended when making frequent requests).
        # Reference: http://cactus.nci.nih.gov/chemical/structure/{molecule name}/smiles
        time.sleep(15)
        return cirpy.resolve(name, 'smiles')
    except Exception as e:
        print(str(e))
        return "FAILED"
