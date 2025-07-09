import os
import logging
import pandas as pd
import csv
import re
from helper import json_extractor, get_output_path, token_formater
from typing import List, Dict, Any, Optional, Callable

# benchmark configuration
from config import (
    MODEL,
    TEMPERATURES,
    BENCHMARKS_ROOT,
    CHALLENGE_FILES,
    FORMULA_OPTIONS,
    PROMPT_TYPES,
    INFERENCE_CALL
)

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

class NMRBenchmarkRunner:
    def __init__(self, benchmarks_root: str, model: str, temperatures: List[float], 
                 challenge_files: List[str], formula_options: List[bool], 
                 prompt_types: Dict[str, Callable], inference_call: Callable):
        self.benchmarks_root = benchmarks_root
        self.model = model
        self.temperatures = temperatures
        self.challenge_files = challenge_files
        self.formula_options = formula_options
        self.prompt_types = prompt_types
        self.inference_call = inference_call
    def clean_model_response(self, response: str) -> str:
        # Clean model response of special tokens and normalize for processing.
        if not response:
            return ""
        cleaned = response.replace('"', '').replace("'", "")
        cleaned = re.sub(r'<[^>]+>', '', cleaned)
        cleaned = re.sub(r'\s+', ' ', cleaned)
        #  special tokens (dependent on model)
        special_tokens = ['<|endoftext|>', '<|startoftext|>', '<pad>', '<unk>', '<s>', '</s>']
        for token in special_tokens:
            cleaned = cleaned.replace(token, '')
        return cleaned.strip().lower()
        
    def run_batch(self, temperature: float, challenge_file: str, outputfile: str, prompt: Callable, formula: bool):
        #Run batch processing for a single configuration guide.
        results = []
        molecules = json_extractor(challenge_file)
        
        for index, molecule in enumerate(molecules):
            logging.info(f"[ID]-{molecule['nmr_challenge_id']}")
            data = token_formater(molecule, formula)
            model_prediction = self.inference_call(
                prompt=prompt(data, formula), 
                temperature=temperature)
            task = {
                "Id": molecule["nmr_challenge_id"],
                "Formula": molecule["formula"],
                "Prediction":self.clean_model_response(model_prediction)
            }
            results.append(task)

        
        # Save raw model response file
        data = pd.DataFrame(results)
        data.to_csv(outputfile, index=False, quoting=csv.QUOTE_ALL)
        logging.info(outputfile)
    
    def run_all_experiments(self) -> None:
        try:
            for prompt_id, prompt_func in self.prompt_types.items():
                # Run entire benchmark
                for formula_option in self.formula_options:
                    for temperature in self.temperatures:
                        for index, challenge_file in enumerate(self.challenge_files):
                            category = (self.challenge_files[index].split("_")[2].split(".json")[0]).upper()
                            file_exist, outputfile = get_output_path(
                                self.benchmarks_root, category, self.model, 
                                prompt_id, temperature, formula_option
                            )
                            
                            logging.info(f"\033[33m\nTest Configuration:\n - Category: {category}\n - Model: {self.model}\n - Prompt: {prompt_id}\n - Temperature: {temperature}\n - Formula: {formula_option}\n \033[0m \n\n")
                            self.run_batch(
                                temperature=temperature,
                                challenge_file=self.challenge_files[index],
                                outputfile=outputfile,
                                prompt=prompt_func,
                                formula=formula_option
                            )
  
        except Exception as e:
            logging.error(f"[ERROR]\n{str(e)}")
            raise

def main() -> None:
    NMRBenchmarkRunner(
        benchmarks_root=BENCHMARKS_ROOT,
        model=MODEL,
        temperatures=TEMPERATURES,
        challenge_files=CHALLENGE_FILES,
        formula_options=FORMULA_OPTIONS,
        prompt_types=PROMPT_TYPES,
        inference_call=INFERENCE_CALL
    ).run_all_experiments()

if __name__ == "__main__":
        main()
    