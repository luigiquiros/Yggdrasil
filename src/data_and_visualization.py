#define function
# Function to process a single CSV file and retrieve Q codes
def process_csv_file(input_file, output_folder):
    # Read the original CSV file containing the species names
    df_species = pd.read_csv(input_file)
    
    # Define the Wikidata Query Service endpoint URL
    endpoint_url = "https://query.wikidata.org/sparql"
    
    # Initialize an empty list to store the Q codes
    q_codes = []
    
    # Initialize a placeholder value for cases where no Q code is found
    placeholder_value = "Not Found"
    
    # Define the number of requests allowed per minute (adjust according to Wikidata's rate limits)
    requests_per_minute = 30
    
    # Iterate through the species names in the CSV file
    for index, row in df_species.iterrows():
        species_name = row[species_header]
    
        # Check if the species_name is a valid string
        if isinstance(species_name, str):
            # Query Wikidata to retrieve the corresponding Q code for the species
            sparql_query = f"""
            SELECT ?species ?speciesLabel WHERE {{
                
                ?species rdfs:label "{species_name}"@en.
                SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }}
            }}
            """
    
            response = requests.get(endpoint_url, headers={'Accept': 'text/csv'}, params={'query': sparql_query})
    
            if response.status_code == 200:
                csv_data = response.text
                df_qcode = pd.read_csv(StringIO(csv_data))
    
                if not df_qcode.empty:
                    q_code = df_qcode['species'].values[0]  # Extract the Q code
                    q_codes.append(q_code)  # Append the Q code to the list
                    print(f"Retrieved Q code {q_code} for species: {species_name}")
                else:
                    q_codes.append(placeholder_value)  # Add a placeholder value
                    print(f"No Q code found for species: {species_name}")
            else:
                q_codes.append(placeholder_value)  # Add a placeholder value for failed requests
                print(f"Failed to retrieve Q code for species: {species_name}. Status code: {response.status_code}")
    
            # Implement rate limiting by waiting between requests
            time.sleep(60 / requests_per_minute)  # Wait for one minute divided by the allowed requests per minute
        else:
            q_codes.append(placeholder_value)  # Add a placeholder value for invalid species names
    
    # Add the Q codes hyperlinks to the original CSV with a new column
    df_species['wikidata_Qcode_hyperlink'] = q_codes
    
    # Add the Q codes only to the original CSV with a new column
    df_species['wikidata_Qcode'] = df_species['wikidata_Qcode_hyperlink'].str.rsplit('/', n=1).str[-1]
    
    # Save the updated CSV file with Q codes in the output folder
    output_csv_filename = os.path.join(output_folder, os.path.basename(input_file))
    df_species.to_csv(output_csv_filename, index=False)
    print(f"Saved updated CSV with Q codes to {output_csv_filename}")
    
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
            'wikidata_Qcode_species': item['species']['value'].split('/')[-1],
            'Species': item['speciesLabel']['value']
        })

    return pd.DataFrame(species)

def recover_LOTUS_data(input_file, output_folder):
    # Load the LOTUSDB CSV
    lotusdb_df = pd.read_csv(LOTUSDB, low_memory=False)
    # Load the CSV file with Q codes
    df_species = pd.read_csv(input_file)

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

        output_filename = os.path.join(output_folder, f"{q_code}.tsv")
        grouped_df.to_csv(output_filename, index=False, sep='\t')

        print(f"Saved grouped data for Q code {q_code} to {output_filename}")
   

def process_LOTUS_data(output_folder, input_folder, LOTUSDB_rc):
    """
    Process LOTUS data and update a CSV file with combined information.

    Parameters:
    output_folder (str): Path to the folder containing .tsv files.
    input_folder (str): Path to the folder containing input CSV files.
    LOTUSDB_rc (str): Path to the LOTUS database CSV file.

    Returns:
    pandas.DataFrame: Updated DataFrame with combined information.
    """
    # Step 1: Initialize variables to store data
    total_compounds = {}
    chemical_classes = {}
    chemical_superclasses = {}
    references = {}
    hyperlinks = {}

    # Step 2: Iterate through .tsv files in the folder
    species_data_folder = os.path.join(output_folder)
    for filename in os.listdir(species_data_folder):
        if filename.endswith(".tsv"):
            qcode = filename.split(".")[0]  # Extract Qcode from the filename

            # Load the .tsv file into a DataFrame
            df_compounds = pd.read_csv(os.path.join(species_data_folder, filename), sep='\t')

            # Step 3: Calculate frequencies of chemical classes for each Qcode (excluding 'not classified')
            filtered_classes = df_compounds[df_compounds['structure_taxonomy_npclassifier_03class'] != 'Not Classified']
            grouped = filtered_classes['structure_taxonomy_npclassifier_03class'].value_counts().reset_index()
            frequencies = grouped.apply(lambda row: f"{row.name} {row['structure_taxonomy_npclassifier_03class']}", axis=1)
            chemical_classes[qcode] = "|".join(frequencies)

            # Step 4: Calculate frequencies of chemical superclasses for each Qcode (excluding 'not classified')
            filtered_sclasses = df_compounds[df_compounds['structure_taxonomy_npclassifier_02superclass'] != 'Not Classified']
            sgrouped = filtered_sclasses['structure_taxonomy_npclassifier_02superclass'].value_counts().reset_index()
            frequencies = sgrouped.apply(lambda row: f"{row.name} {row['structure_taxonomy_npclassifier_02superclass']}", axis=1)
            chemical_superclasses[qcode] = "|".join(frequencies)

    # Step 7: Load the general info .csv into a DataFrame
    csv_file = None
    for filename in os.listdir(input_folder):
        if filename.endswith(".csv"):
            csv_file = os.path.join(input_folder, filename)
            break  # Stop searching after finding the first CSV file

    if csv_file is not None:
        # Load the CSV file into the df_general_info DataFrame
        df_general_info = pd.read_csv(csv_file)

        # Step 7: Create new columns for combined information
        df_general_info['predicted_class'] = df_general_info['wikidata_Qcode'].map(chemical_classes)
        df_general_info['predicted_superclass'] = df_general_info['wikidata_Qcode'].map(chemical_superclasses)

        # Step 7: Save the updated DataFrame to the .csv file
        output_csv_file = os.path.join(input_folder, 'Full_results.csv')
        df_general_info.to_csv(output_csv_file, index=False, sep=',')

        # Replace values matching the pattern with an empty cell
        df_general_info['predicted_class'] = df_general_info['predicted_class'].replace('index|predicted_class', '')
        df_general_info['predicted_superclass'] = df_general_info['predicted_superclass'].replace('index|predicted_superclass', '')

        # Add reported compounds to each Q code
        LotusDB_rc = pd.read_csv(LOTUSDB_rc)
        spDB = LotusDB_rc[['wikidata_Qcode', 'Reported_comp_Species']].drop_duplicates()
        df = pd.merge(df_general_info, spDB, how='left', left_on='wikidata_Qcode', right_on='wikidata_Qcode')
        df.drop('wikidata_Qcode', axis=1, inplace=True)
        # Save the final DataFrame
        df.to_csv(output_csv_file, index=False, sep=',')
        print("CSV file updated in the input folder.")
    
    else:
        print("No CSV file found in the input folder.")
    
    return df


# Define base colors for each Pathway
pathway_shades= {
    'Terpenoids': ('#618264', '#D0E7D2'),  # Green start and lighter green end
    'Alkaloids': ('#305F72', '#5CBCE2'),   # Blue start and lighter blue end  
    'Shikimates and Phenylpropanoids': ('#80558C', '#CBA0AE'),  # Purple start and lighter purple end
    'Polyketides': ('#EF4B4B', '#EC8F6A'),  # Red start and lighter purple end
    'Fatty acids': ('#FF6C22', '#FF9209'),  # Orange start and lighter purple end
    'Amino acids and Peptides': ('#F4E869', '#FAF2D3'),  # Yellow start and lighter purple end
    'Carbohydrates': ('#65451F','#C8AE7D')  # Brown start and lighter purple end
}

# Define custom colors for the 7 pathway categories
pathway_colors = {
    'Terpenoids': '#618264',  # Green start and lighter green end
    'Alkaloids': '#305F72',   # Blue start and lighter blue end
    'Shikimates and Phenylpropanoids': '#80558C',  # Purple start and lighter purple end
    'Polyketides': '#EF4B4B',  # Red start and lighter purple end
    'Fatty acids': '#FF6C22',  # Orange start and lighter purple end
    'Amino acids and Peptides': '#F4E869',  # Yellow start and lighter purple end
    'Carbohydrates': '#65451F' # Brown start and lighter purple end
}

def interpolate_color(color1, color2, factor: float):
    """Interpolate between two colors"""
    color1 = np.array(mcolors.to_rgb(color1))
    color2 = np.array(mcolors.to_rgb(color2))
    return mcolors.to_hex((1 - factor) * color1 + factor * color2)

def generate_shades(pathway, num_shades):
    base_color, end_color = pathway_shades.get(pathway, ('gray', 'lightgray'))
    if num_shades == 1:
        return [base_color]  # Return the base color if only one shade is requested
    shades = []
    for i in range(num_shades):
        factor = i / (num_shades - 1)
        shades.append(interpolate_color(base_color, end_color, factor))
    return shades

def split_chemical_superclass(row):
    # Check if the value is a string before splitting
    if isinstance(row['chemical_superclass'], str):
        parts = row['chemical_superclass'].split('-')
        if len(parts) == 2:
            return parts[0], parts[1]  # Pathway and Superclass are present
        else:
            return parts[0], 'Unknown'  # Only Pathway is present, or the format is not as expected
    else:
        # Return default values if the entry is NaN or not a string
        return 'Unknown', 'Unknown'

def barplot_sclass_species(output_folder):
    """
    Process data from .tsv files in the output folder and visualize it with a stacked barplot.

    Parameters:
    output_folder (str): Path to the folder containing .tsv files.

    Returns:
    None
    """
    # Step 1: Read data from all .tsv files in the output_folder
    all_data = pd.concat([pd.read_csv(os.path.join(output_folder, filename), sep='\t') 
                          for filename in os.listdir(output_folder) if filename.endswith(".tsv")])

    # Step 2: Rename the "organism_taxonomy_09species" column to "species"
    all_data.rename(columns={'organism_taxonomy_09species': 'species'}, inplace=True)

    # Remove 'API Error-API Error' and 'Not Classified-Not Classified'
    all_data = all_data[~all_data['chemical_superclass'].isin(['API Error-API Error', 'Not Classified-Not Classified'])]

    # Step 3: Apply the function to each row
    all_data[['Pathway', 'superclass']] = all_data.apply(lambda row: split_chemical_superclass(row), axis=1, result_type='expand')

    # Step 4: Process data for color mapping
    color_map = {}
    for pathway, superclasses in all_data.groupby('Pathway')['superclass'].unique().items():
        shades = generate_shades(pathway, len(superclasses))
        for superclass, shade in zip(superclasses, shades):
            color_map[f"{pathway}-{superclass}"] = shade

    # Step 5: Group and aggregate data to calculate recurrence
    agg_data = all_data.groupby(['species', 'chemical_superclass']).size().reset_index(name='recurrence')

    # Convert 'species' column to categorical data
    agg_data['species'] = pd.Categorical(agg_data['species'], categories=agg_data['species'].unique(), ordered=True)

    # Get unique species names
    unique_species = agg_data['species'].unique()

    # Get unique chemical superclasses and sort them alphabetically
    unique_superclasses = sorted(agg_data['chemical_superclass'].unique())

    # Calculate total recurrence for each species
    total_recurrence_per_species = agg_data.groupby('species')['recurrence'].sum()

    # Create the stacked barplot with the custom color palette
    fig = px.bar(agg_data, y='species', x='recurrence',
                 title='Stacked Barplot of Predicted Superclasses Occurrence for Species',
                 labels={'recurrence': 'Recurrence'},
                 color='chemical_superclass',
                 color_discrete_map=color_map,
                 category_orders={'species': unique_species, 'chemical_superclass': unique_superclasses},
                 orientation='h')

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Species<i>')

    # Set species labels in italics
    fig.update_layout(yaxis=dict(tickmode='array',
                                  tickvals=list(range(len(unique_species))),
                                  ticktext=[f'<i>{species}</i>' for species in unique_species]
                                  ))

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Modify the size of the figure
    fig.update_layout(width=1500, height=1500)

    # Save the figure as an HTML file
    fig.write_html(f'{output_folder}Wikidata_superclass_barplot_species.html')

    # Show the figure
    fig.show()


def barplot_sclass_species_normalized(output_folder):
    """
    Process data from .tsv files in the output folder and visualize it with a normalized stacked barplot.

    Parameters:
    output_folder (str): Path to the folder containing .tsv files.

    Returns:
    None
    """
    # Step 1: Read data from all .tsv files in the output_folder
    all_data = pd.concat([pd.read_csv(os.path.join(output_folder, filename), sep='\t') 
                          for filename in os.listdir(output_folder) if filename.endswith(".tsv")])

    # Step 2: Rename the "organism_taxonomy_09species" column to "species"
    all_data.rename(columns={'organism_taxonomy_09species': 'species'}, inplace=True)

    # Remove 'API Error-API Error' and 'Not Classified-Not Classified'
    all_data = all_data[~all_data['chemical_superclass'].isin(['API Error-API Error', 'Not Classified-Not Classified'])]

    # Apply the function to each row
    all_data[['Pathway', 'superclass']] = all_data.apply(lambda row: split_chemical_superclass(row), axis=1, result_type='expand')

    # Step 4: Process data for color mapping
    color_map = {}
    for pathway, superclasses in all_data.groupby('Pathway')['superclass'].unique().items():
        shades = generate_shades(pathway, len(superclasses))
        for superclass, shade in zip(superclasses, shades):
            color_map[f"{pathway}-{superclass}"] = shade

    # Step 5: Group and aggregate data to calculate recurrence
    agg_data = all_data.groupby(['species', 'chemical_superclass']).size().reset_index(name='recurrence')

    # Normalize the recurrence values within each species group
    agg_data['recurrence_normalized'] = agg_data.groupby('species')['recurrence'].transform(lambda x: x / x.sum()) * 100

    # Convert 'species' column to categorical data
    agg_data['species'] = pd.Categorical(agg_data['species'], categories=agg_data['species'].unique(), ordered=True)

    # Get unique species names
    unique_species = agg_data['species'].unique()

    # Get unique chemical superclasses and sort them alphabetically
    unique_superclasses = sorted(agg_data['chemical_superclass'].unique())

    # Calculate total recurrence for each species
    total_recurrence_per_species = agg_data.groupby('species')['recurrence'].sum()

    # Create the stacked barplot with the custom color palette
    fig = px.bar(agg_data, y='species', x='recurrence_normalized',
                 title='Normalized Stacked Barplot of Predicted Superclasses Occurrence for Species',
                 labels={'recurrence_normalized': 'Recurrence'},
                 color='chemical_superclass',
                 color_discrete_map=color_map,
                 category_orders={'species': unique_species, 'chemical_superclass': unique_superclasses},
                 orientation='h')

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Species<i>')

    # Set species labels in italics
    fig.update_layout(yaxis=dict(tickmode='array',
                                  tickvals=list(range(len(unique_species))),
                                  ticktext=[f'<i>{species}</i>' for species in unique_species]
                                  ))

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Get the maximum value of the x-axis
    max_x = agg_data['recurrence_normalized'].max()

    for species, total_recurrence in total_recurrence_per_species.items():
        # Find the maximum recurrence for this species
        max_recurrence = agg_data[agg_data['species'] == species]['recurrence_normalized'].max()

        # Position the annotation at the end of the bar with a slight offset
        x_position = max_x + 2
        if x_position > 100:  # Ensure the label is not outside the plot area
            x_position = 100  # Set it to 100 if it exceeds the maximum x-axis value

        # Add annotation to the plot
        fig.add_annotation(
            x=x_position,
            y=species,
            text=f'Total compounds: {total_recurrence}',
            showarrow=False,
            font=dict(size=10, color='black'),
            xanchor='left',
            yanchor='middle'
        )

    # Modify the size of the figure
    fig.update_layout(width=1500, height=1500)

    # Save the figure as an HTML file
    fig.write_html(f'{output_folder}Wikidata_superclass_barplot_species_normalized.html')

    # Show the figure
    fig.show()
