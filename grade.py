import os
import re
import json
import pandas as pd
import csv
import logging 
import subprocess
from typing import List, Dict, Any, Optional
from smile import generate_smiles_string 
from tanimoto_similarity import calculate_tanimoto, canonical_smiles_match

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

from config import (
    MODEL,
    BENCHMARKS_ROOT,
    CHALLENGE_IDS,
    GRADED_DIR
)

class NMRGrader:
    def __init__(self, model: str, benchmarks_root: str, challenge_ids: List[str], graded_dir: str):
        self.model = model
        self.benchmarks_root = benchmarks_root
        self.challenge_ids = challenge_ids
        self.graded_dir = graded_dir
        self.max_text_length = 20000
        
    def csv_extractor_local(self, filename: str) -> List[Dict[str, Any]]:
        try:
            file = pd.read_csv(filename)
            molecules = pd.DataFrame(file)
            
            return [{
                'Id': Id,
                "Formula": Formula,                              
                "Prediction": Prediction,           
            } 
            for Id, Formula, Prediction
            in zip(molecules['Id'], molecules['Formula'], molecules['Prediction'])
            ]
        except Exception as e:
            logging.error(f"Error extracting CSV {filename}: {str(e)}")
            raise
    
    def get_paths(self, root: str) -> Optional[List[str]]:
        """Get all paths in a directory."""
        try:
            com = subprocess.run(["ls", f"{root}"], stdout=subprocess.PIPE, text=True).stdout
            return None if com == ""  else [f"{root}/{x}" for x in com.strip().split("\n")]    
        except Exception as e:
            logging.error(f"Paths issues: {e}")
            raise

    def truncate_prediction(self, text: str) -> str:
        if len(text) <= self.max_text_length:
            return text
        keep_length = (self.max_text_length - 3) // 2
        return text[:keep_length] + "..." + text[-keep_length:]


    def extract_answer(self, text: str) -> str:
        """Extract molecule name from structured model response."""
        if not text:
            raise ValueError("model response is empty check inference call")
        match = re.search(r'###\s*Start answer\s*###\s*(.*?)\s*###\s*End answer\s*###', text, re.IGNORECASE | re.DOTALL)
        if not match:
            return "FAILED EXTRACTION"
        raw = match.group(1).strip()
        cleaned = re.sub(r'[^a-zA-Z0-9\s\-\+\(\),]', '', raw)
        cleaned = re.sub(r'\s+', ' ', cleaned).strip()
        return cleaned if cleaned else "FAILED EXTRACTION"

    def verdict_calculation(self, name_match: bool, smiles_correct: str, smiles_llm: str, tanimoto: float) -> int:
        smiles_match = canonical_smiles_match(smiles_correct, smiles_llm)
        exact_tanimoto = (tanimoto == 1.0)
        return 1 if (name_match or smiles_match or exact_tanimoto) else 0


    def grade_task_file(self, task_file: str, reference_hash: List[Dict[str, Any]]) -> None:
        graded_tasks = []
        tasks = self.csv_extractor_local(task_file)
        for task in tasks:
            if task['Id'] in reference_hash:
                entry = reference_hash[task['Id']].copy()
                entry["Prediction"] = self.extract_answer(task['Prediction'])   
                name_match = entry["TrueName"].lower().strip() == entry["Prediction"].lower().strip()
                entry["SMILE-LLM"] = entry['SMILE-correct'] if name_match else generate_smiles_string(entry["Prediction"])                
                entry["TanimotoCoefficient"] = calculate_tanimoto(
                    entry['SMILE-correct'], 
                    entry["SMILE-LLM"]
                )
                entry["Verdict"] = self.verdict_calculation(
                    name_match,
                    entry['SMILE-correct'],
                    entry["SMILE-LLM"],
                    entry["TanimotoCoefficient"]
                )
                entry["ScratchPad"] = self.truncate_prediction(task['Prediction'])
                graded_tasks.append(entry)
        
        # Save graded results
        data = pd.DataFrame(graded_tasks)
        outputfile = task_file.replace(str(self.benchmarks_root), str(self.graded_dir))
        os.makedirs(os.path.dirname(outputfile), exist_ok=True)
        data.to_csv(outputfile, index=False, quoting=csv.QUOTE_ALL)
        logging.info(outputfile)

    def load_correct_references(self) -> Dict[str, Dict[str, Any]]:
        """Load all reference data."""
        refs = {}
        for challenge_id in self.challenge_ids:
            df = pd.read_csv(challenge_id)
            category = challenge_id.split('_')[-1].replace(".csv", "").upper()
            refs[category] = {
                row['NMR_Challenge_ID']: {
                    'Id': row['NMR_Challenge_ID'],
                    'Formula': row['Formula'],
                    'TrueName': row['True names'],
                    'SMILE-correct': row['Smiles']
                }
                for _, row in df.iterrows()
            }
        return refs

    def get_raw_files(self) -> List[Dict[str, Any]]:
        """Get all raw files to grade organized by category."""
        raw_files = []
        
        model_path = f"{self.benchmarks_root}/{self.model}"
        categories = self.get_paths(model_path)
        for category in categories:
            entry = {
                'category': category.split('/')[-1],
                'task_files': []
            }
            formula_paths = self.get_paths(category)
            if formula_paths:
                for formula_path in formula_paths:
                    task_paths = self.get_paths(formula_path)
                    if task_paths:
                        entry['task_files'].extend(task_paths)
            raw_files.append(entry)
        return raw_files

    def grade_all_files(self) -> None:
        """Grade all files in the benchmark directory."""
        try:
            # Load correct references
            correct_references = self.load_correct_references()
            # Get files to grade
            raw_files = self.get_raw_files()
            for category_data in raw_files:
                category = category_data['category']
                if category not in correct_references:
                    logging.warning(f"No correct references found for category: {category}")
                    continue
                
                reference_hash = correct_references[category]
                for task_file in category_data['task_files']:
                    self.grade_task_file(task_file, reference_hash)
                    
        except Exception as e:
            logging.error(f"Error in grade_all_files: {str(e)}")
            raise

def main() -> None:
    try:
        grader = NMRGrader(
            model=MODEL,
            benchmarks_root=BENCHMARKS_ROOT,
            challenge_ids=CHALLENGE_IDS,
            graded_dir=GRADED_DIR
        )
        
        grader.grade_all_files()
        
    except Exception as e:
        logging.error(f"Fatal error in main: {str(e)}")
        raise

if __name__ == "__main__":
    main()