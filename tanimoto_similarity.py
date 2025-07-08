from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, DataStructs
from rdkit.Chem import rdMolDescriptors

def calculate_tanimoto(smile1, smile2):
    try:
        mol1 = Chem.MolFromSmiles(smile1)
        mol2 = Chem.MolFromSmiles(smile2)

        if mol1 is None or mol2 is None:
            print(f"[ERROR-1] Invalid SMILES: '{smile1}' or '{smile2}'")
            return 0

        fpg = rdFingerprintGenerator.GetMorganGenerator()
        fp1 = fpg.GetFingerprint(mol1)
        fp2 = fpg.GetFingerprint(mol2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    
    except Exception as e:
        print(f"[ERROR-0] Tanimoto coefficient calculation failed: {e}")
        return 0
