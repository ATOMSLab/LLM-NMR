# Benchmark Evaluation Instructions

## Overview

The prompts in this folder were used for comparison against human performance evaluation.

## Data

This benchmark run uses the **raw, unnormalized full H + CNMR data**.

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
