system_prompt="""You are an expert in NMR spectroscopy analysis and organic chemistry."""

def base_prompt(nmr_data, formula=None):
    """Zero-shot prompt"""
    return f"""
            Here is 1H NMR spectrum data for a certain molecule {":" if formula is None else f"with chemical formula {formula}:"}
            {nmr_data}
            What's the molecule's name?
            {"" if formula is None else "Be mindful of stoichiometry and ensure consistency with the given formula." }
            Format the final answer like this - 
            ### Scratchpad ### <scratchpad> ### Scratchpad ###
            ### Start answer ### <prediction> ### End answer ###
            The prediction should only contain the name of the molecule and no other text
        """

def cot_prompt(nmr_data, formula=None):
    """Chain of thought prompt"""
    return f"""
            Here is 1H NMR spectrum data for a certain molecule {":" if formula is None else f"with chemical formula {formula}:"}
            {nmr_data}
            What's the molecule's name?
            Think step-by-step, making extensive use of a scratchpad to record your thoughts.
            Consider finding ways to group related peaks together, and keep track of the stoichiometry 
            and the amount of unassigned H atoms as you make provisional assignments.
            {"" if formula is None else "Be mindful of stoichiometry and ensure consistency with the given formula." }
            Format the final answer like this - 
            ### Scratchpad ### <scratchpad> ### Scratchpad ###
            ### Start answer ### <prediction> ### End answer ###
            The prediction should only contain the name of the molecule and no other text
        """

def logic_tips_prompt(nmr_data, formula=None):
    """Chain of thought prompt with logic tips"""
    return f"""
           Here is 1H NMR spectrum data for a certain molecule {":" if formula is None else f"with chemical formula {formula}:"}
            {nmr_data}
            What's the molecule's name?
            Let's think step-by-step. 
            Use a scratchpad to record your thoughts. First, identify all groups of peaks and the corresponding number of hydrogens, and ensure the stoichiometry of the formula matches that in the spectrum Then, for each group of peaks, hypothesize one or more molecule fragments that could explain the signals in that peak.
            Rank the hypotheses according to which is most plausible
            Occasionally synthesize a candidate molecule from the provisionally-assigned fragments
            Be mindful of connectivity and ensure consistency with all hypothesized fragments
            {"" if formula is None else "Be mindful of stoichiometry and ensure consistency with the given formula." }
            Continue using the scratchpad until you're confident about the answer.
            Format the final answer like this - 
            ### Scratchpad ### <scratchpad> ### Scratchpad ###
            ### Start answer ### <prediction> ### End answer ###
            The prediction should only contain the name of the molecule and no other text
        """

def expert_tips_prompt(nmr_data, formula=None):
    """Chain of thought prompt with expert tips"""
    return f"""
            Here is 1H NMR spectrum data for a certain molecule {":" if formula is None else f"with chemical formula {formula}:"}
            {nmr_data} 
            What's the molecule's name? Let's think step-by-step. 
            Use a scratchpad to record your thoughts, and continuing using it until your confident about your answer. 
            
            Before the analysis of NMR spectra, think about the molecular formula. Is the molecule fully saturated, or are there some multiple bonds or rings?
            Consider the solvent used. Do you expect the appearance of signals from exchangeable hydrogens (e.g., in CDCl3or DMSO), or are they exchanged to deuterium (e.g. in D2O or methanol-d4)? 
            Identify solvent signals and ignore them in the structural analysis.
            Consider the number of signals, their positions (chemical shifts), intensity, shape and splitting.
            The number of signals corresponds to the number of non-equivalent atoms in the molecule.
            Typical regions in proton NMR spectra are: carboxylic acids: 10–13 ppm, aldehydes: 9–11 ppm, aromatic hydrogens: 6.5–8 ppm, hydrogens attached to a double-bond carbon: 4–7.5 ppm, O–CH fragments (alcohols, ethers): 3–5 ppm, hydrogens attached to a triple-bond carbon: 2–3.5 ppm, hydrogens attached to a single-bond carbon without an electronegative substituent: 0.5–2.5 ppm.
            
            {"" if formula is None else "Be mindful of stoichiometry and ensure consistency with the given formula." }  
            Format the final answer like this - 
            ### Scratchpad ### <scratchpad> ### Scratchpad ###
            ### Start answer ### <prediction> ### End answer ###
            The prediction should only contain the name of the molecule and no other text
        """

def expert_logic_tips_prompt(nmr_data, formula=None):
    """Chain of thought prompt with both expert and logic tips"""
    return f"""
            Here is 1H NMR spectrum data for a certain molecule {":" if formula is None else f"with chemical formula {formula}:"}
            {nmr_data} 
            Let's think step-by-step using both logical analysis and expert knowledge.
            
            Use a scratchpad to:
            1. Identify all peak groups in the 1H NMR spectrum and assign them to their corresponding hydrogens.
            2. Ensure stoichiometry matches the formula
            3. For each peak group, hypothesize possible molecular fragments
            4. Rank hypotheses by plausibility
            5. Build candidate molecules from fragments
            
            Consider these expert guidelines:
            - Check if the molecule is fully saturated or has multiple bonds/rings
            - Look for exchangeable hydrogens based on solvent
            - Note typical 1H NMR chemical shift regions:
              * 10-13 ppm: carboxylic acids
              * 9-11 ppm: aldehydes
              * 6.5-8 ppm: aromatic hydrogens
              * 4-7.5 ppm: alkene hydrogens
              * 3-5 ppm: O-CH fragments
              * 2-3.5 ppm: alkyne hydrogens
              * 0.5-2.5 ppm: alkyl hydrogens
            
            {"" if formula is None else "Be mindful of stoichiometry and ensure consistency with the given formula." }
            Format the final answer like this - 
            ### Scratchpad ### <scratchpad> ### Scratchpad ###
            ### Start answer ### <prediction> ### End answer ###
            The prediction should only contain the name of the molecule and no other text
        """

# self augmentation prompt
def consistency_reprompt():
    return """
            is the predicted molecule name consistent with the formula, answer with 
            ### Scratchpad ###
            ### Start answer ###<yes or no>### End answer ###
            """


def regeneration_prompt():
    return f"""
            re-analyse the spectra data
            ### Scratchpad ### <scratchpad> ### Scratchpad ###
            ### Start answer ### <prediction> ### End answer ###
            The prediction should only contain the name of the molecule and no other text
            """
