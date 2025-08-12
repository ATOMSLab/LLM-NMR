# Benchmark Evaluation Instructions

## Overview

The prompts in this folder were used for comparison against human performance evaluation.

## Data

This benchmark run uses the raw, unnormalized full H + CNMR data, which can be found in the unzipped `.datasets/human_comparison_benchmark folder.`
It includes data for both CNMR and HNMR.
You will need to provide a function that sorts the tasks and aggregates the two data modalities before prompt injection.

## Steps to Recreate

1. **Prepare Prompts**
   Create a new file named `prompts_CNMR.py` and add the prompt cases based on the `.txt` files in this directory.

1. **Update Configuration**
   Edit `config.py` to set the correct dataset path,prompts file and output path.

1. **Install Dependencies**
   Install the required packages using the following command:

   ```bash
   pip install -r requirements.txt
   ```

1. **Run the Benchmark**
   Execute the main script:

   ```bash
   python main.py
   ```
