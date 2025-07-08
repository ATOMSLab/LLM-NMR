import os
import logging
from pathlib import Path
from helper import json_extractor, csv_extractor, write_benchmark_result, get_output_path, token_formater
from tanimoto_similarity import calculate_tanimoto
from smile import generate_smiles_string

from config import (
    MODEL,
    TEMPERATURES,
    BENCHMARKS_ROOT,
    CHALLENGE_FILES,
    CHALLENGE_IDS,
    FORMULA_OPTIONS,
    PROMPT_TYPES,
    INFERENCE_CALL

)

def grade():
    pass
       

def run_batch(temperature,id_file,challenge_file,outputfile,prompt,self_aug,formula):
    correct_reference=csv_extractor(id_file)
    molecules= json_extractor(challenge_file)
    for index,molecule in enumerate (molecules):
        logging.info(f"[ID]-{molecule['nmr_challenge_id']}")
        data=token_formater(molecule, formula)
        formula = molecule.get('formula') if formula else None
        response,model_prediction=INFERENCE_CALL(prompt=prompt(data,formula),temperature=temperature)
        scratch_pad=response.split("### Start answer ###")[0].replace("\n"," ")
        #grade call 




if (__name__=="__main__"):
 try:
    for ID,PROMPT in PROMPT_TYPES.items():
    #run entire benchmark
        for FORMULA in FORMULA_OPTIONS:
            for TEMP in TEMPERATURES:
                for index,file in enumerate(CHALLENGE_FILES):
                    CATEGORY=(CHALLENGE_FILES[index].split("_")[2].split(".json")[0]).upper()
                    verdict,OUTPUTFILE=get_output_path(BENCHMARKS_ROOT,CATEGORY, MODEL, ID, TEMP, FORMULA)
                    logging.info(f"\033[33m\nTest Configuration:\n - Category: {CATEGORY}\n - Model: {MODEL}\n - Prompt: {ID}\n - Temperature: {TEMP}\n - Formula: {FORMULA}\n \033[0m \n\n")
                    run_batch(
                            temperature=TEMP,
                            id_file=CHALLENGE_IDS[index],
                            challenge_file=CHALLENGE_FILES[index],
                            outputfile=OUTPUTFILE,
                            prompt=PROMPT,
                            formula=FORMULA
                        )

 except Exception as e:
        logging.error("[ERROR]\n",str(e))
        raise
        