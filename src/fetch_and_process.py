import os
import pandas as pd
import requests

def fetch_species_from_qcode(qcode):
    """
    Fetches all species under a given genus using the Wikidata SPARQL endpoint.

    Parameters:
    genus_qcode (str): The Wikidata Q-code for the genus.

    Returns:
    pd.DataFrame: A DataFrame containing the species and their corresponding Q-codes.
    """

    # SPARQL query to fetch species under a given genus
    query = """
    SELECT ?species ?speciesLabel WHERE {
      ?species wdt:P171* wd:%s .
      ?species wdt:P105 wd:Q7432.
      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
    }
    """ % qcode

    # URL for the Wikidata SPARQL endpoint
    url = "https://query.wikidata.org/sparql"

    # Request headers
    headers = {
        "User-Agent": "Wikidata Species Fetcher/0.1 (https://www.wikidata.org/wiki/Wikidata:Data_access)"
    }

    # Perform the request
    response = requests.get(url, headers=headers, params={'query': query, 'format': 'json'})

    if response.status_code != 200:
        raise Exception("Failed to fetch data: HTTP status code {}".format(response.status_code))

    # Parse the JSON response
    data = response.json()

    # Extract results
    species = []
    for item in data['results']['bindings']:
        species.append({
            'wikidata_Qcode': item['species']['value'].split('/')[-1],
            'Species': item['speciesLabel']['value']
        })

    return pd.DataFrame(species)

def recover_LOTUS_data_sp(input_file, lotusdb_path, output_folder):
    # Load the LOTUSDB CSV
    lotusdb_df = pd.read_csv(lotusdb_path, low_memory=False)
    # Load the CSV file with Q codes
    df_species = pd.read_csv(input_file)

    # Create a subfolder for species data
    sp_subfolder = os.path.join(output_folder, 'species_data')
    os.makedirs(sp_subfolder, exist_ok=True)

    # Iterate through the Q codes in the species CSV
    for q_code in df_species['wikidata_Qcode']:  # Ensure this column name matches your species CSV
        # Skip processing if Q code is 'Not Found'
        if q_code == 'Not Found':
            continue

        # Filter LOTUSDB data for the current Q code
        filtered_lotusdb = lotusdb_df[lotusdb_df['wikidata_Qcode'] == q_code]  # Ensure this column name matches your LOTUSDB CSV

        # Group and aggregate the data
        grouped_df = filtered_lotusdb.groupby("structure_inchikey").agg({
            # Add all the aggregation rules here
            # Example:
            "structure_wikidata": "first",
            "structure_inchi": "first",
            "structure_smiles": "first",
            "structure_molecular_formula": "first",
            "structure_exact_mass": "first",
            "structure_xlogp": "first",
            "structure_smiles_2D": "first",
            "structure_cid": "first",
            "structure_nameIupac": "first",
            "structure_nameTraditional": "first",
            "structure_taxonomy_npclassifier_01pathway": "first",
            "structure_taxonomy_npclassifier_02superclass": "first",
            "structure_taxonomy_npclassifier_03class": "first",
            "organism_wikidata": "first",
            "organism_taxonomy_gbifid": "first",
            "organism_taxonomy_ncbiid": "first",
            "organism_taxonomy_ottid": "first",
            "organism_taxonomy_01domain": "first",
            "organism_taxonomy_02kingdom": "first",
            "organism_taxonomy_03phylum": "first",
            "organism_taxonomy_04class": "first",
            "organism_taxonomy_05order": "first",
            "organism_taxonomy_06family": "first",
            "organism_taxonomy_07tribe": "first",
            "organism_taxonomy_08genus": "first",
            "organism_taxonomy_09species": "first",
            "organism_taxonomy_10varietas": "first",
            "reference_wikidata": lambda x: "|".join(map(str, x)),
            "reference_doi": lambda x: "|".join(map(str, x))
        }).reset_index()
        
        # Create 'chemical_superclass' and 'chemical_class' columns
        grouped_df['chemical_superclass'] = grouped_df['structure_taxonomy_npclassifier_01pathway'] + '-' + grouped_df['structure_taxonomy_npclassifier_02superclass']
        grouped_df['chemical_class'] = grouped_df['structure_taxonomy_npclassifier_01pathway'] + '-' + grouped_df['structure_taxonomy_npclassifier_03class']
        
        # Save the grouped data as a TSV file with the Q code as the filename
        #output_subfolder = os.path.join(output_folder, 'species_data')
        #os.makedirs(output_subfolder, exist_ok=True)

        output_filename = os.path.join(sp_subfolder, f"{q_code}.tsv")
        grouped_df.to_csv(output_filename, index=False, sep='\t')

        print(f"Saved grouped data for Q code {q_code} to {output_filename}")
   
def recover_LOTUS_data_g(input_file, lotusdb_path, output_folder):
    # Load the LOTUSDB CSV
    lotusdb_df = pd.read_csv(lotusdb_path, low_memory=False)
    # Load the CSV file with Q codes
    df_species = pd.read_csv(input_file)

    # Create a subfolder for genus data
    genus_subfolder = os.path.join(output_folder, 'genus_data')
    os.makedirs(genus_subfolder, exist_ok=True)

    # Extract unique genus values from the 'Genus' column
    unique_genus = df_species['Genus'].unique()

    # Iterate through unique genus values
    for genus in unique_genus:
        # Skip processing if Q code is 'Not Found'
        if genus == 'Not Found':
            continue

        # Filter LOTUSDB data for the current Q code
        filtered_lotusdb = lotusdb_df[lotusdb_df['organism_taxonomy_08genus'] == genus]  # Ensure this column name matches your LOTUSDB CSV

        # Group and aggregate the data
        grouped_df = filtered_lotusdb.groupby("structure_inchikey").agg({
            # Add all the aggregation rules here
            # Example:
            "structure_wikidata": "first",
            "structure_inchi": "first",
            "structure_smiles": "first",
            "structure_molecular_formula": "first",
            "structure_exact_mass": "first",
            "structure_xlogp": "first",
            "structure_smiles_2D": "first",
            "structure_cid": "first",
            "structure_nameIupac": "first",
            "structure_nameTraditional": "first",
            "structure_taxonomy_npclassifier_01pathway": "first",
            "structure_taxonomy_npclassifier_02superclass": "first",
            "structure_taxonomy_npclassifier_03class": "first",
            "organism_wikidata": "first",
            "organism_taxonomy_gbifid": "first",
            "organism_taxonomy_ncbiid": "first",
            "organism_taxonomy_ottid": "first",
            "organism_taxonomy_01domain": "first",
            "organism_taxonomy_02kingdom": "first",
            "organism_taxonomy_03phylum": "first",
            "organism_taxonomy_04class": "first",
            "organism_taxonomy_05order": "first",
            "organism_taxonomy_06family": "first",
            "organism_taxonomy_07tribe": "first",
            "organism_taxonomy_08genus": "first",
            "organism_taxonomy_09species": "first",
            "organism_taxonomy_10varietas": "first",
            "reference_wikidata": lambda x: "|".join(map(str, x)),
            "reference_doi": lambda x: "|".join(map(str, x))
        }).reset_index()
        
        # Create 'chemical_superclass' and 'chemical_class' columns
        grouped_df['chemical_superclass'] = grouped_df['structure_taxonomy_npclassifier_01pathway'] + '-' + grouped_df['structure_taxonomy_npclassifier_02superclass']
        grouped_df['chemical_class'] = grouped_df['structure_taxonomy_npclassifier_01pathway'] + '-' + grouped_df['structure_taxonomy_npclassifier_03class']
        
        # Save the grouped data as a TSV file in the genus_data subfolder
        output_filename = os.path.join(genus_subfolder, f"{genus}.tsv")
        grouped_df.to_csv(output_filename, index=False, sep='\t')

        print(f"Saved grouped data for genus {genus} to {output_filename}")


def process_species_data(input_folder, output_folder, lotusdb_path):
    """
    Processes species-level data from .tsv files, calculates frequency of chemical classes and superclasses,
    merges with general species information, and adds reported compound counts.

    Parameters:
    - input_folder: Path to the folder containing species general info CSV file.
    - output_folder: Path where processed results will be saved.
    - LOTUSDB_path: Path to the LOTUSDB CSV file.
    """

    # Step 1: Initialize dictionaries to store data
    chemical_classes = {}
    chemical_superclasses = {}

    # Step 2: Iterate through .tsv files in the folder
    species_data_folder = os.path.join(output_folder, 'species_data')
    if not os.path.exists(species_data_folder):
        print(f"Error: {species_data_folder} does not exist.")
        return

    for filename in os.listdir(species_data_folder):
        if filename.endswith(".tsv"):
            qcode = filename.split(".")[0]  # Extract Qcode from the filename

            # Load the .tsv file into a DataFrame
            df_compounds = pd.read_csv(os.path.join(species_data_folder, filename), sep='\t')

            # Step 3: Calculate frequencies of chemical classes
            if 'chemical_class' in df_compounds.columns:
                class_counts = df_compounds['chemical_class'].value_counts()
                chemical_classes[qcode] = "|".join([f"{count} {cls}" for cls, count in class_counts.items()])

            # Step 4: Calculate frequencies of chemical superclasses
            if 'chemical_superclass' in df_compounds.columns:
                superclass_counts = df_compounds['chemical_superclass'].value_counts()
                chemical_superclasses[qcode] = "|".join([f"{count} {scls}" for scls, count in superclass_counts.items()])

    # Step 5: Load the general info CSV file
    csv_file = None
    for filename in os.listdir(input_folder):
        if filename.endswith(".csv"):
            csv_file = os.path.join(input_folder, filename)
            break  # Stop searching after finding the first CSV file

    if csv_file is None:
        print("Error: No general info CSV file found in the input folder.")
        return

    df_general_info = pd.read_csv(csv_file)

    # Step 6: Create new columns for combined information
    df_general_info['predicted_class'] = df_general_info['wikidata_Qcode'].map(chemical_classes)
    df_general_info['predicted_superclass'] = df_general_info['wikidata_Qcode'].map(chemical_superclasses)

    # Step 7: Load the LOTUSDB CSV
    LOTUSDB = pd.read_csv(lotusdb_path, low_memory=False)

    # Step 8: Merge LOTUSDB data to add reported compounds for each Q code
    df = pd.merge(df_general_info, LOTUSDB[['wikidata_Qcode', 'Reported_comp_Species', 'Reported_comp_Genus' ]],
                  how='left', left_on='wikidata_Qcode', right_on='wikidata_Qcode')

    # Step 9: Save the updated DataFrame to a CSV file
    output_csv_file = os.path.join(output_folder, 'Full_results.csv')
    df.to_csv(output_csv_file, index=False, sep=',')

    print(f"✅ Process completed! Results saved to {output_csv_file}")

