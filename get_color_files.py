import pandas as pd
import os
import glob
# Please note the this script ONLY work if the phage-genome-analysis-pipeline snakefile was used before
# Please enter to path to the result, This path should be the same as the analysis_folder_path in phage-genome-analysis-pipeline snakefile 
gbk_path = ["/home/champa/DATA/PHAGE_genome_analysis/Champoux_A/Analysis_05-2026/Analysis/"]
# Create a list to stock all the dataframe
dataframes = []

# Iterating throught all the "phold_per_cds_predictions.tsv" file. There should be one file per sample 
# (one of those file is created in each phold analysis)
for dir in gbk_path:
    for file_path in glob.glob(os.path.join(dir, "*/Pharokka/*cds_final_merged_output.tsv")):
        print(file_path)

        # Thean each of the tsv file is converted into a pandas dataframe
        df = pd.read_csv(file_path, delimiter='\t')
        
        # Only the needed columns are selected
        selected_columns = df[['gene','category']] 

        # QC: print(selected_columns)
        # Those columns are added to the list previoulsy created
        dataframes.append(selected_columns)

# Concatenating all the dataframe (2 columns only)
final_df = pd.concat(dataframes)

# the dataframe is exported in a csv file
final_df.to_csv("merged_cds_predictions.csv", index=False, sep=',')

# A hexadecimal code of a colour is attribute to every function found by phold 
function_to_color = {
    'head and packaging': '#900c3f',
    'unknown function': '#fb4fd9',
    'connector':'#8f2be7',
    'tail': '#007bd8',
    'lysis': '#00e1da',
    'other':'#1fb819',
    'DNA, RNA and nucleotide metabolism': '#ffdb28',
    'moron, auxiliary metabolic gene and host takeover': '#f28200',
    'transcription regulation': '#03fc9d',
    'integration and excision':'#0d5c0a'
}

# Extracting unique functions and their associated colors
unique_functions = final_df['category'].unique()
print(unique_functions)
#If a function doesn't have an attribute colour, the arow will be black
color_mapping = [(func, function_to_color.get(func, '#000000')) for func in unique_functions]

# Create a DataFrame for the function-to-color mapping
color_df = pd.DataFrame(color_mapping, columns=['function', 'color'])

# Exporting a csv file with the fonction and the colours
color_df.to_csv("function_color_mapping.csv", index=False, sep=',')