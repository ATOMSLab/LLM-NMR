
# LLM-NMR: Benchmarking LLMs for NMR Spectral Reasoning

[![arXiv](https://img.shields.io/badge/arXiv-2025.XXXXX-b31b1b.svg)](https://arxiv.org/abs/XXXXX)

**LLM-NMR** is a benchmark framework for evaluating large language models (LLMs) on NMR spectral analysis tasks, assessing their chemical reasoning capabilities using 1D NMR data.

---

## 🔍 Overview

This repo includes:

* ✅ An extensive dataset of 115 NMR problems with varying difficulty levels
* 🧠 An Inference script for running LLMs from public facing apis for OpenAI, Anthropic, and Google
* 📊 Tools for grading model outputs using SMILES + Tanimoto comparison
* 🔁 Experimental configuration for prompts, temperature, formula masking, etc.

---

## 📦 Installation

```bash
git clone https://github.com/ATOMSLab/LLM-NMR.git
cd LLM-NMR
pip install -r requirements.txt
```

**Requirements:**
`pandas`, `numpy`, `rdkit`, `cirpy`, `openai`, `anthropic`, `google-generativeai`, `json`, `matplotlib`, `seaborn`

---

## ⚙️ Setup

### 1. Configure `.env`

Add your API keys:

```dotenv
OPENAI_API_KEY=...
ANTHROPIC_API_KEY=...
GOOGLE_API_KEY=...
```

### 2. Edit `config.py`

Set:

* Path to dataset JSON
* Output directory for inference and grading
* Model settings (name, temperature, etc.)

---

## 🚀 Run Benchmark

### Step 1: Run Inference

```bash
python3 main.py
```

Generates output CSV (configurable path), e.g.:

```csv
Id,Formula,Prediction
191,C4H8O2,"### Scratchpad ###...### Start answer ###Ethyl acetate### End answer ###"
```

### Step 2: Grade Outputs

```bash
python3 grade.py
```

This will:

* Extract final answer from LLM output
* Convert ground truth and prediction to SMILES
* Compute Tanimoto similarity
* Aggregate results for score

---

## 📁 Dataset

* Download JSON dataset from [Drive link]()
* Place it in the root directory (or update path in `config.py`)

### Composition:

| Difficulty | Count |
| ---------- | ----- |
| Easy       | 53    |
| Medium     | 38    |
| Hard       | 24    |



---

## 🧪 Models Tested

| Model             | Provider  | Type      |
| ----------------- | --------- | --------- |
| GPT-4o            | OpenAI    | Standard  |
| GPT-4o mini       | OpenAI    | Standard  |
| o1                | OpenAI    | Reasoning |
| o1-mini           | OpenAI    | Reasoning |
| o3-mini           | OpenAI    | Reasoning |
| Claude-3.5 Sonnet | Anthropic | Standard  |
| Gemini-2.0-Flash  | Google    | Standard  |

---

## 🔍 Prompting Modes

| Strategy | Description                    |
| -------- | ------------------------------ |
| P1       | Minimal instruction            |
| P2       | Chain-of-Thought (CoT)         |
| P3       | CoT + domain logic             |
| P4       | CoT + expert NMR tips          |
| P5       | CoT + knowledge + logic (full) |

---

## ⚗️ Experiment Variables

* **Temperature:** `0.0`, `0.5`, `0.8`, `1.0`
* **Formula Inclusion:** with/without `molecular_formula`
* **Reasoning Effort:** Low / Medium / High (for o-models)
* **Difficulty Tier:** Easy / Medium / Hard

---

## 📈 Results (With Formula)

| Model             | Accuracy | Rank |
| ----------------- | -------- | ---- |
| o1                | 69%      | 🥇   |
| o3-mini           | 65%      | 🥈   |
| Claude-3.5 Sonnet | 51%      | 🥉   |
| Gemini-2.0-Flash  | 38%      | 4th  |
| GPT-4o            | 25%      | 5th  |
| o1-mini           | 30%      | 6th  |
| GPT-4o mini       | 10%      | 7th  |

---

## 📚 Citation

```bibtex
@article{llm_spectroscopy_2025,
  title={LLM Spectroscopy: A Benchmark for Chemical Reasoning with AI},
  author={[Authors]},
  journal={...},
  year={2025}
}
```

---

## 🔗 Links

* 🔬 [NMR-Challenge.com](https://nmr-challenge.com) – source of benchmark problems
* 📂 [Download Dataset]() – (link to JSON dataset on Google Drive)
* 🧪 [RDKit Documentation](https://www.rdkit.org/docs/) – for SMILES handling and similarity metrics


