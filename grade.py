
import os
import time
import logging
import warnings
import pandas as pd
from rdkit import Chem, __version__ as rdkit_version
from rdkit.Chem import DataStructs, rdFingerprintGenerator



class NMRGrader:
    """NMR Grader focused on grading functionality"""
    
    def __init__(self):
        self.logger = setup_logging()
        self.logger.info("NMR Grader Initialized")
        

        self.processed_files = []
        self.error_files = []
        

        self.answer_keys = {}
        self._load_answer_keys()
    
    
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
    
    def calculate_verdict(self, row):
        # Name match
        if 'prediction' in row and 'true_name' in row:
            pred_name = str(row['prediction']).lower().strip() if pd.notna(row['prediction']) else ""
            true_name = str(row['true_name']).lower().strip() if pd.notna(row['true_name']) else ""
            
            if pred_name and true_name and pred_name == true_name:
                return 1
    
        if 'llm_smiles' in row and 'true_smiles' in row:
            if pd.notna(row['llm_smiles']) and pd.notna(row['true_smiles']):
                if str(row['llm_smiles']) == str(row['true_smiles']):
                    return 1
        
        if pd.notna(row.get('CanonicalSMILES_Match')) and row['CanonicalSMILES_Match']:
            return 1
        
        if pd.notna(row.get('TanimotoCoefficient')) and row['TanimotoCoefficient'] > 0.99:
            if 'llm_smiles' in row and pd.notna(row['llm_smiles']):
                llm_smiles = str(row['llm_smiles'])
                if '.' in llm_smiles:
                    return 1
                else:
                    if pd.notna(row.get('CanonicalSMILES_Match')):
                        return 1 if row['CanonicalSMILES_Match'] else 0
                    else:
                        return 0
            else:
                return 1
        
        return 0
    
    def grade_file(self, file_path, difficulty=None):
        
        try:
            df = pd.read_csv(file_path)
            original_columns = df.columns.tolist()
            self.logger.info(f"Grading file: {file_path} ({len(df)} rows)")
            
            self.logger.debug(f"Original columns: {original_columns}")
            
            columns_removed = []
            for col in GRADING_COLUMNS_TO_REMOVE:
                if col in df.columns:
                    df = df.drop(columns=[col])
                    columns_removed.append(col)
            
            
            
            df['id'] = df['id'].astype(str)
            
            answer_df = self.answer_keys[difficulty].copy()
            if 'id' in answer_df.columns:
                answer_df['id'] = answer_df['id'].astype(str)
            
            answer_cols_needed = ['id']
            if 'true_name' in answer_df.columns:
                answer_cols_needed.append('true_name')
            if 'true_smiles' in answer_df.columns:
                answer_cols_needed.append('true_smiles')
            
            answer_df_subset = answer_df[answer_cols_needed]
            
            merged_df = pd.merge(df, answer_df_subset, 
                               on='id', how='left', suffixes=('', '_answer'))
            
            print(f"  Calculating Tanimoto coefficients...")
            merged_df['TanimotoCoefficient'] = merged_df.apply(
                lambda row: self.calculate_tanimoto(row['llm_smiles'], row['true_smiles']),
                axis=1
            )
            
            merged_df['CanonicalSMILES_Match'] = merged_df.apply(
                lambda row: self.canonical_smiles_match(row['llm_smiles'], row['true_smiles']),
                axis=1
            )
            
            merged_df['Verdict'] = merged_df.apply(self.calculate_verdict, axis=1)
            
            total = len(merged_df)
            correct = merged_df['Verdict'].sum()
            accuracy = (correct / total * 100) if total > 0 else 0
            
            reverse_mapping = {}
            for orig_col in original_columns:
                for standard, variations in COLUMN_MAPPINGS.items():
                    if orig_col in variations:
                        if standard in merged_df.columns:
                            reverse_mapping[standard] = orig_col
                        break
            
            merged_df = merged_df.rename(columns=reverse_mapping)
            
            final_columns = []
            for col in original_columns:
                if col in merged_df.columns and col not in GRADING_COLUMNS_TO_REMOVE:
                    final_columns.append(col)
            
            grading_cols_to_add = []
            
            if 'True names' in original_columns or any('True names' in col for col in original_columns):
                grading_cols_to_add.append('True names')
                if 'true_name' in merged_df.columns:
                    merged_df = merged_df.rename(columns={'true_name': 'True names'})
            elif 'TrueName' in original_columns:
                grading_cols_to_add.append('TrueName')
                if 'true_name' in merged_df.columns:
                    merged_df = merged_df.rename(columns={'true_name': 'TrueName'})
            elif 'true_name' in merged_df.columns:
                grading_cols_to_add.append('true_name')
            
            if 'Smiles' in original_columns or 'SMILES' in original_columns:
                col_name = 'SMILES' if 'SMILES' in original_columns else 'Smiles'
                grading_cols_to_add.append(col_name)
                if 'true_smiles' in merged_df.columns:
                    merged_df = merged_df.rename(columns={'true_smiles': col_name})
            elif 'SMILE-correct' in original_columns:
                grading_cols_to_add.append('SMILE-correct')
                if 'true_smiles' in merged_df.columns:
                    merged_df = merged_df.rename(columns={'true_smiles': 'SMILE-correct'})
            elif 'true_smiles' in merged_df.columns:
                grading_cols_to_add.append('true_smiles')
            
            grading_cols_to_add.extend(['TanimotoCoefficient', 'CanonicalSMILES_Match', 'Verdict'])
            
            for col in grading_cols_to_add:
                if col in merged_df.columns and col not in final_columns:
                    final_columns.append(col)
            
            merged_df = merged_df[final_columns]
            
            merged_df.to_csv(file_path, index=False)
            self.logger.info(f"Saved graded file: {file_path}")
            
            self.processed_files.append({
                'path': file_path,
                'difficulty': difficulty,
                'total': total,
                'correct': correct,
                'accuracy': accuracy
            })
            
            print(f"  ✓ Graded successfully - Accuracy: {accuracy:.1f}% ({correct}/{total})")
            
            return merged_df
            
        except Exception as e:
            self.logger.error(f"Error grading file {file_path}: {e}")
            self.error_files.append({'path': file_path, 'error': str(e)})
            print(f"  Error: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def grade_directory_standard(self):
        """Grade files using standard directory structure"""
        try:
            model_dirs = [d for d in os.listdir(BENCHMARK_DIR) 
                         if os.path.isdir(os.path.join(BENCHMARK_DIR, d))]
            
            if not model_dirs:
                print(f"No model directories found in {BENCHMARK_DIR}")
                return False
            
            print(f"\nFound {len(model_dirs)} model directories")
            
            for model_idx, model_name in enumerate(model_dirs):
                model_dir = os.path.join(BENCHMARK_DIR, model_name)
                print(f"\n=== Processing Model {model_idx+1}/{len(model_dirs)}: {model_name} ===")
                
                for difficulty in ["EASY", "MEDIUM", "HARD"]:
                    if difficulty not in self.answer_keys:
                        continue
                    
                    for formula_condition in FORMULA_CONDITIONS:
                        dir_path = os.path.join(model_dir, difficulty, formula_condition)
                        
                        if not os.path.exists(dir_path):
                            continue
                        
                        csv_files = [f for f in os.listdir(dir_path) if f.endswith('.csv')]
                        
                        for csv_file in csv_files:
                            file_path = os.path.join(dir_path, csv_file)
                            print(f"\nGrading: {model_name}/{difficulty}/{formula_condition}/{csv_file}")
                            self.grade_file(file_path, difficulty)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error in standard grading: {e}")
            return False
# In[ ]:





# This will overwrite the CSV files! 

# In[ ]:


def run_grading():
    """Main function to run grading"""
    print("\n" + "="*60)
    print("NMR GRADING SYSTEM")
    print("="*60)
    print(f"\nBenchmark directory: {BENCHMARK_DIR}")
    print(f"Answer key directory: {ANSWER_KEY_DIR}")
    print(f"Formula conditions: {FORMULA_CONDITIONS}")
    
    grader = NMRGrader()
    
    start_time = time.time()
    
    success = grader.grade_directory_standard()
    
    elapsed = time.time() - start_time
    
    grader.print_summary()
    
    print(f"\nTotal time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    print("\nGrading complete!")
    
    return grader 

print("Main execution function defined!")


# In[ ]:





# In[ ]:



grader = run_grading()

