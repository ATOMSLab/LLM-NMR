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

def canonical_smiles_match(self, smiles1, smiles2):
        """Check if two SMILES represent the same molecule"""
        if not rdkit_available:
            return False
        
        if not smiles1 or not smiles2 or pd.isna(smiles1) or pd.isna(smiles2):
            return False
        
        try:
            mol1 = Chem.MolFromSmiles(str(smiles1))
            mol2 = Chem.MolFromSmiles(str(smiles2))
            
            if mol1 is None or mol2 is None:
                return False
            
            can1 = Chem.MolToSmiles(mol1, canonical=True, isomericSmiles=True)
            can2 = Chem.MolToSmiles(mol2, canonical=True, isomericSmiles=True)
            
            return can1 == can2
            
        except Exception as e:
            self.logger.debug(f"Error comparing canonical SMILES: {e}")
            return False