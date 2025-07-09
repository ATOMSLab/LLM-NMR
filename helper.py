import json
import pandas as pd
import os

def json_extractor(filename):
 try:
    file = open(filename,"r")
    txt=file.read()
    json_data=json.loads(txt)
    return json_data
 except Exception as e:
    print("[ERROR]\n",str(e))
    
def csv_extractor(filename):
    try:
        file=pd.read_csv(filename)
        molecules=pd.DataFrame(file)
        return [({
            'id':id,
            'formula':formula.strip(),
            'name':name.lower().strip(),
            'smile':smile
        }) 
                for id,formula,name,smile in zip(molecules['NMR_Challenge_ID'],molecules['Formula'],molecules['True names'],molecules['Smiles'])]
    except Exception as e:
        print("[ERROR]\n",str(e))


def token_formater(data, include_formula):
    try:
        peaks = []
        for peak in data["peaks"]:
            peak_str = "--------\n"
            peak_str += f"Name: {peak['name']}\n"
            peak_str += f"Shift: {peak['shift']}\n"
            peak_str += f"Range: {' '.join(str(value) for value in peak['range'])}\n"
            if include_formula:  # Only include hydrogens if formula is included
                peak_str += f"Hydrogens: {peak['hydrogens']}\n"
            peak_str += f"Integral: {peak['integral']}\n"
            peak_str += f"Class: {peak['class']}\n"
            peak_str += f"J values: {' '.join(str(value) for value in peak['j_values'])}\n"
            
            if 'method' in peak:
                peak_str += f"Method: {peak['method']}\n"
                
            peaks.append(peak_str)
        
        if include_formula:
            formatted_data = f"\nFormula: {data.get('formula', 'N/A')}\nPeaks:\n{''.join(peaks)}"
        else:
            formatted_data = f"\nPeaks:\n{''.join(peaks)}"

        return formatted_data
        
    except Exception as e:
        logging.error(f"Error formatting data: {str(e)}")
        return ""
    

def prepare_output_file(path):
    try:
        if (os.path.exists(path)):
            return
        else :
            csv_file=open(path,"a")
            csv_file.write("Id,Formula,TrueName,SMILE-correct,Prediction,SMILE-LLM,Verdict,TanimotoCoefficient,Isomers,ScratchPad")
    except Exception as e:
        print(str(e))

# Generate output file path based on parameters
def get_output_path(root, difficulty, model, prompt_type, temperature, with_formula):
    ensure_benchmark_dirs(root, model)
    form_type = "WITH_FORMULA" if with_formula else "NO_FORMULA"
    base_path = root / model / difficulty / form_type
    filename = f"{difficulty}-results-T{temperature}-{prompt_type}-{model}.csv"
    path=base_path / filename
    if os.path.exists(path):
        return True,path
    prepare_output_file(path)
    return False,path

# create benchmark directories if they don't exist
def ensure_benchmark_dirs(root, model):
    for difficulty in ["EASY", "MEDIUM", "HARD"]:
        for form_type in ["WITH_FORMULA", "NO_FORMULA"]:
            path = root/ model / difficulty / form_type
            path.mkdir(parents=True, exist_ok=True)
