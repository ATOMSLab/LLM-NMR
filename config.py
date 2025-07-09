# benchmark config
import os
from pathlib import Path
from prompts import base_prompt, cot_prompt, logic_tips_prompt, expert_tips_prompt, expert_logic_tips_prompt
from llms import call_openAI, call_gemini, call_claude_sonnet


# Model configuration
MODEL = "gemini-2.0-flash-exp"
INFERENCE_CALL=call_gemini
TEMPERATURES = [0, 0.5, 0.8, 1.0]

# Paths
IDS_ROOT = "./H1 NMR Datasets/IDS/"
CHALLENGE_ROOT = "./H1 NMR Datasets/Challenges/"
BENCHMARKS_ROOT = Path("./BENCHMARK")

#files 
CHALLENGE_FILES = [
    f"{CHALLENGE_ROOT}Hnmr_spectra_easy.json",
    f"{CHALLENGE_ROOT}Hnmr_spectra_medium.json",
    f"{CHALLENGE_ROOT}Hnmr_spectra_hard.json",
]

CHALLENGE_IDS = [
    f"{IDS_ROOT}Hnmr_challenge_easy_ids.csv",
    f"{IDS_ROOT}Hnmr_challenge_medium_ids.csv",
    f"{IDS_ROOT}Hnmr_challenge_hard_ids.csv",
]

# Formula options
FORMULA_OPTIONS = [True, False]
# prompts
PROMPT_TYPES = {
    "base": base_prompt,
    "cot": cot_prompt,
    "logic": logic_tips_prompt,
    "expert": expert_tips_prompt,
    "expert_logic": expert_logic_tips_prompt
}
