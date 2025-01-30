import os
import numpy as np
import pandas as pd
import plotly.express as px
import matplotlib.colors as mcolors

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
 
 
def plot_species_superclass(output_folder, split_chemical_superclass, generate_shades):
    """
    Reads .tsv files, processes species and superclass data, and generates a stacked bar plot.

    Parameters:
    - output_folder: Path to the folder containing 'species_data' subfolder.
    - split_chemical_superclass: Function to split 'chemical_superclass' into 'Pathway' and 'Superclass'.
    - generate_shades: Function to generate color shades for pathways.
    
    Saves the plot as an HTML file in the output folder.
    """

    species_data_folder = os.path.join(output_folder, 'species_data')

    # Step 1: Read data from all .tsv files
    all_data = pd.concat([
        pd.read_csv(os.path.join(species_data_folder, filename), sep='\t')
        for filename in os.listdir(species_data_folder) if filename.endswith(".tsv")
    ], ignore_index=True)

    # Step 2: Rename 'organism_taxonomy_09species' to 'species'
    all_data.rename(columns={'organism_taxonomy_09species': 'species'}, inplace=True)

    # Step 3: Apply the function to split 'chemical_superclass'
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

    # Get unique species names and superclasses
    unique_species = agg_data['species'].unique()
    unique_superclasses = sorted(agg_data['chemical_superclass'].unique())

    # Step 6: Create the stacked barplot
    fig = px.bar(
        agg_data, y='species', x='recurrence',
        title='Stacked Barplot of Predicted Superclasses Occurrence for Species',
        labels={'recurrence': 'Recurrence'},
        color='chemical_superclass',
        color_discrete_map=color_map,
        category_orders={'species': unique_species, 'chemical_superclass': unique_superclasses},
        orientation='h'
    )

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Species<i>')

    # Set species labels in italics
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=list(range(len(unique_species))),
        ticktext=[f'<i>{species}</i>' for species in unique_species]
    ))

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Modify the figure size
    fig.update_layout(width=1500, height=1500)

    # Save the figure as an HTML file
    output_html_file = os.path.join(output_folder, 'Wikidata_superclass_barplot_species.html')
    fig.write_html(output_html_file)

    # Show the figure
    fig.show()

    print(f"✅ Process completed! Visualization saved to {output_html_file}")



def plot_species_superclass_norm(output_folder, split_chemical_superclass, generate_shades):
    """
    Reads .tsv files, processes species and superclass data, normalizes recurrence values, and generates a stacked bar plot.

    Parameters:
    - output_folder: Path to the folder containing 'species_data' subfolder.
    - split_chemical_superclass: Function to split 'chemical_superclass' into 'Pathway' and 'Superclass'.
    - generate_shades: Function to generate color shades for pathways.

    Saves the plot as an HTML file in the output folder.
    """

    species_data_folder = os.path.join(output_folder, 'species_data')

    # Step 1: Read data from all .tsv files
    all_data = pd.concat([
        pd.read_csv(os.path.join(species_data_folder, filename), sep='\t')
        for filename in os.listdir(species_data_folder) if filename.endswith(".tsv")
    ], ignore_index=True)

    # Step 2: Rename 'organism_taxonomy_09species' to 'species'
    all_data.rename(columns={'organism_taxonomy_09species': 'species'}, inplace=True)

    # Step 3: Remove specific values
    all_data = all_data[~all_data['chemical_superclass'].isin(['API Error-API Error', 'Not Classified-Not Classified'])]

    # Step 4: Apply function to split 'chemical_superclass'
    all_data[['Pathway', 'superclass']] = all_data.apply(lambda row: split_chemical_superclass(row), axis=1, result_type='expand')

    # Step 5: Process data for color mapping
    color_map = {}
    for pathway, superclasses in all_data.groupby('Pathway')['superclass'].unique().items():
        shades = generate_shades(pathway, len(superclasses))
        for superclass, shade in zip(superclasses, shades):
            color_map[f"{pathway}-{superclass}"] = shade

    # Step 6: Group and aggregate data to calculate recurrence
    agg_data = all_data.groupby(['species', 'chemical_superclass']).size().reset_index(name='recurrence')

    # Step 7: Normalize the recurrence values within each species group
    agg_data['recurrence_normalized'] = agg_data.groupby('species')['recurrence'].transform(lambda x: x / x.sum()) * 100

    # Convert 'species' column to categorical data
    agg_data['species'] = pd.Categorical(agg_data['species'], categories=agg_data['species'].unique(), ordered=True)

    # Get unique species names and superclasses
    unique_species = agg_data['species'].unique()
    unique_superclasses = sorted(agg_data['chemical_superclass'].unique())

    # Step 8: Calculate total recurrence for each species
    total_recurrence_per_species = agg_data.groupby('species')['recurrence'].sum()

    # Step 9: Create the stacked barplot
    fig = px.bar(
        agg_data, y='species', x='recurrence_normalized',
        title='Normalized Stacked Barplot of Predicted Superclasses Occurrence for Species',
        labels={'recurrence_normalized': 'Normalized Recurrence (%)'},
        color='chemical_superclass',
        color_discrete_map=color_map,
        category_orders={'species': unique_species, 'chemical_superclass': unique_superclasses},
        orientation='h'
    )

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Species<i>')

    # Set species labels in italics
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=list(range(len(unique_species))),
        ticktext=[f'<i>{species}</i>' for species in unique_species]
    ))

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Step 10: Add annotations for total compounds per species
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

    # Step 11: Modify figure size
    fig.update_layout(width=1500, height=1500)

    # Save the figure as an HTML file
    output_html_file = os.path.join(output_folder, 'Wikidata_superclass_barplot_species_normalized.html')
    fig.write_html(output_html_file)

    # Show the figure
    fig.show()

    print(f"✅ Process completed! Visualization saved to {output_html_file}")


def plot_species_pathway(output_folder):
    """
    Process data from .tsv files in the output folder and visualize it with a stacked barplot.

    Parameters:
    - output_folder (str): Path to the folder containing 'species_data' subfolder.

    Saves the visualization as an HTML file.
    """

    # Define custom colors for the 7 pathway categories
    pathway_colors = {
        'Terpenoids': '#618264',  # Green
        'Alkaloids': '#305F72',   # Blue
        'Shikimates and Phenylpropanoids': '#80558C',  # Purple
        'Polyketides': '#EF4B4B',  # Red
        'Fatty acids': '#FF6C22',  # Orange
        'Amino acids and Peptides': '#F4E869',  # Yellow
        'Carbohydrates': '#65451F'  # Brown
    }

    species_data_folder = os.path.join(output_folder, 'species_data')

    # Step 1: Read data from all .tsv files
    all_data = pd.concat([
        pd.read_csv(os.path.join(species_data_folder, filename), sep='\t')
        for filename in os.listdir(species_data_folder) if filename.endswith(".tsv")
    ], ignore_index=True)

    # Step 2: Rename columns
    all_data.rename(columns={'organism_taxonomy_09species': 'species'}, inplace=True)
    all_data.rename(columns={'structure_taxonomy_npclassifier_01pathway': 'Pathway'}, inplace=True)

    # Step 3: Remove 'API Error' and 'Not Classified' from 'Pathway'
    all_data = all_data[~all_data['Pathway'].isin(['API Error', 'Not Classified'])]

    # Step 4: Group and aggregate data to calculate recurrence
    agg_data = all_data.groupby(['species', 'Pathway']).size().reset_index(name='recurrence')

    # Convert 'species' column to categorical data
    agg_data['species'] = pd.Categorical(agg_data['species'], categories=agg_data['species'].unique(), ordered=True)

    # Get unique species names and pathways
    unique_species = agg_data['species'].unique()
    unique_pathways = sorted(agg_data['Pathway'].unique())

    # Step 5: Create the stacked barplot with custom colors
    fig = px.bar(
        agg_data, y='species', x='recurrence',
        title='Stacked Barplot of Predicted Pathways Occurrence for Species',
        labels={'recurrence': 'Recurrence'},
        color='Pathway',
        color_discrete_map=pathway_colors,  # Use custom colors
        category_orders={'species': unique_species, 'Pathway': unique_pathways},
        orientation='h'
    )

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Species<i>')

    # Set species labels in italics
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=list(range(len(unique_species))),
        ticktext=[f'<i>{species}</i>' for species in unique_species]
    ))

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Modify the figure size
    fig.update_layout(width=1500, height=1500)

    # Save the figure as an HTML file
    output_html_path = os.path.join(output_folder, 'Wikidata_pathway_barplot_species.html')
    fig.write_html(output_html_path)

    # Show the figure
    fig.show()

    print(f"✅ Process completed! Visualization saved to {output_html_path}")

def plot_species_pathway_norm(output_folder):
    """
    Reads .tsv files, processes species and pathway data, normalizes recurrence values, and generates a stacked bar plot.

    Parameters:
    - output_folder (str): Path to the folder containing 'species_data' subfolder.

    Saves the plot as an HTML file in the output folder.
    """

    # Define custom colors for the 7 pathway categories
    pathway_colors = {
        'Terpenoids': '#618264',  
        'Alkaloids': '#305F72',  
        'Shikimates and Phenylpropanoids': '#80558C',  
        'Polyketides': '#EF4B4B',  
        'Fatty acids': '#FF6C22',  
        'Amino acids and Peptides': '#F4E869',  
        'Carbohydrates': '#65451F'  
    } 

    species_data_folder = os.path.join(output_folder, 'species_data')

    # Step 1: Read data from all .tsv files
    all_data = pd.concat([
        pd.read_csv(os.path.join(species_data_folder, filename), sep='\t')
        for filename in os.listdir(species_data_folder) if filename.endswith(".tsv")
    ], ignore_index=True)

    # Step 2: Rename columns
    all_data.rename(columns={'organism_taxonomy_09species': 'species'}, inplace=True)
    all_data.rename(columns={'structure_taxonomy_npclassifier_01pathway': 'Pathway'}, inplace=True)

    # Step 3: Remove 'API Error' and 'Not Classified' from 'Pathway'
    all_data = all_data[~all_data['Pathway'].isin(['API Error', 'Not Classified'])]

    # Step 4: Group and aggregate data to calculate recurrence
    agg_data = all_data.groupby(['species', 'Pathway']).size().reset_index(name='recurrence')

    # Step 5: Normalize the recurrence values within each species group
    agg_data['recurrence_normalized'] = agg_data.groupby('species')['recurrence'].transform(lambda x: x / x.sum()) * 100

    # Convert 'species' column to categorical data
    agg_data['species'] = pd.Categorical(agg_data['species'], categories=agg_data['species'].unique(), ordered=True)

    # Get unique species names and pathways
    unique_species = agg_data['species'].unique()
    unique_pathways = sorted(agg_data['Pathway'].unique())

    # Step 6: Calculate total recurrence for each species
    total_recurrence_per_species = agg_data.groupby('species')['recurrence'].sum()

    # Step 7: Create the stacked barplot with custom colors
    fig = px.bar(
        agg_data, y='species', x='recurrence_normalized',
        title='Normalized Stacked Barplot of Predicted Pathways Occurrence for Species',
        labels={'recurrence_normalized': 'Normalized Recurrence (%)'},
        color='Pathway',
        color_discrete_map=pathway_colors,  # Use custom colors
        category_orders={'species': unique_species, 'Pathway': unique_pathways},
        orientation='h'
    )

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Species<i>')

    # Set species labels in italics
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=list(range(len(unique_species))),
        ticktext=[f'<i>{species}</i>' for species in unique_species]
    ))

    # Step 8: Add annotations for total compounds per species
    max_x = agg_data['recurrence_normalized'].max()

    for species, total_recurrence in total_recurrence_per_species.items():
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

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Modify the figure size
    fig.update_layout(width=1500, height=1500)

    # Save the figure as an HTML file
    output_html_path = os.path.join(output_folder, 'Wikidata_pathway_barplot_normalized.html')
    fig.write_html(output_html_path)

    # Show the figure
    fig.show()

    print(f"✅ Process completed! Visualization saved to {output_html_path}")


def plot_genus_superclass(output_folder, generate_shades):
    """
    Reads .tsv files, processes genus and superclass data, and generates a stacked bar plot.

    Parameters:
    - output_folder (str): Path to the folder containing 'genus_data' subfolder.
    - generate_shades: Function to generate color shades for pathways.

    Saves the plot as an HTML file in the output folder.
    """

    genus_data_folder = os.path.join(output_folder, 'genus_data')

    # Step 1: Read data from all .tsv files
    all_data = pd.concat([
        pd.read_csv(os.path.join(genus_data_folder, filename), sep='\t')
        for filename in os.listdir(genus_data_folder) if filename.endswith(".tsv")
    ], ignore_index=True)

    # Step 2: Rename columns
    all_data.rename(columns={'organism_taxonomy_08genus': 'genus'}, inplace=True)

    # Step 3: Remove 'API Error' and 'Not Classified' from 'chemical_superclass'
    all_data = all_data[~all_data['chemical_superclass'].isin(['API Error-API Error', 'Not Classified-Not Classified'])]

    # Function to split 'chemical_superclass' into pathway and superclass
    def split_chemical_superclass(row):
        if isinstance(row['chemical_superclass'], str):
            parts = row['chemical_superclass'].split('-')
            return parts if len(parts) == 2 else [parts[0], 'Unknown']
        return ['Unknown', 'Unknown']

    # Step 4: Apply the function to split chemical superclass
    all_data[['Pathway', 'superclass']] = all_data.apply(lambda row: split_chemical_superclass(row), axis=1, result_type='expand')

    # Step 5: Process data for color mapping
    color_map = {}
    for pathway, superclasses in all_data.groupby('Pathway')['superclass'].unique().items():
        shades = generate_shades(pathway, len(superclasses))
        for superclass, shade in zip(superclasses, shades):
            color_map[f"{pathway}-{superclass}"] = shade

    # Step 6: Group and aggregate data to calculate recurrence
    agg_data = all_data.groupby(['genus', 'chemical_superclass']).size().reset_index(name='recurrence')

    # Convert 'genus' column to categorical data
    agg_data['genus'] = pd.Categorical(agg_data['genus'], categories=agg_data['genus'].unique(), ordered=True)

    # Get unique genus names and chemical superclasses
    unique_genera = agg_data['genus'].unique()
    unique_superclasses = sorted(agg_data['chemical_superclass'].unique())

    # Step 7: Create the stacked barplot with the custom color palette
    fig = px.bar(
        agg_data, y='genus', x='recurrence',
        title='Stacked Barplot of Predicted Superclasses Occurrence for Genus',
        labels={'recurrence': 'Recurrence'},
        color='chemical_superclass',
        color_discrete_map=color_map,  # Use custom colors
        category_orders={'genus': unique_genera, 'chemical_superclass': unique_superclasses},
        orientation='h'
    )

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Genus<i>')

    # Set genus labels in italics
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=list(range(len(unique_genera))),
        ticktext=[f'<i>{genus}</i>' for genus in unique_genera]
    ))

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Modify the figure size
    fig.update_layout(width=1500, height=1600)

    # Save the figure as an HTML file
    output_html_path = os.path.join(output_folder, 'Wikidata_superclass_barplot_genus.html')
    fig.write_html(output_html_path)

    # Show the figure
    fig.show()

    print(f"✅ Process completed! Visualization saved to {output_html_path}")

def plot_genus_superclass_norm(output_folder, generate_shades):
    """
    Reads .tsv files, processes genus and superclass data, normalizes recurrence values, and generates a stacked bar plot.

    Parameters:
    - output_folder (str): Path to the folder containing 'genus_data' subfolder.
    - generate_shades: Function to generate color shades for pathways.

    Saves the plot as an HTML file in the output folder.
    """

    genus_data_folder = os.path.join(output_folder, 'genus_data')

    # Step 1: Read data from all .tsv files
    all_data = pd.concat([
        pd.read_csv(os.path.join(genus_data_folder, filename), sep='\t')
        for filename in os.listdir(genus_data_folder) if filename.endswith(".tsv")
    ], ignore_index=True)

    # Step 2: Rename columns
    all_data.rename(columns={'organism_taxonomy_08genus': 'genus'}, inplace=True)

    # Step 3: Remove 'API Error' and 'Not Classified' from 'chemical_superclass'
    all_data = all_data[~all_data['chemical_superclass'].isin(['API Error-API Error', 'Not Classified-Not Classified'])]

    # Function to split 'chemical_superclass' into pathway and superclass
    def split_chemical_superclass(row):
        if isinstance(row['chemical_superclass'], str):
            parts = row['chemical_superclass'].split('-')
            return parts if len(parts) == 2 else [parts[0], 'Unknown']
        return ['Unknown', 'Unknown']

    # Step 4: Apply the function to split chemical superclass
    all_data[['Pathway', 'superclass']] = all_data.apply(lambda row: split_chemical_superclass(row), axis=1, result_type='expand')

    # Step 5: Process data for color mapping
    color_map = {}
    for pathway, superclasses in all_data.groupby('Pathway')['superclass'].unique().items():
        shades = generate_shades(pathway, len(superclasses))
        for superclass, shade in zip(superclasses, shades):
            color_map[f"{pathway}-{superclass}"] = shade

    # Step 6: Group and aggregate data to calculate recurrence
    agg_data = all_data.groupby(['genus', 'chemical_superclass']).size().reset_index(name='recurrence')

    # Step 7: Normalize the recurrence values within each genus group
    agg_data['recurrence_normalized'] = agg_data.groupby('genus')['recurrence'].transform(lambda x: x / x.sum()) * 100

    # Convert 'genus' column to categorical data
    agg_data['genus'] = pd.Categorical(agg_data['genus'], categories=agg_data['genus'].unique(), ordered=True)

    # Get unique genus names and chemical superclasses
    unique_genera = agg_data['genus'].unique()
    unique_superclasses = sorted(agg_data['chemical_superclass'].unique())

    # Step 8: Calculate total recurrence for each genus
    total_recurrence_per_genus = agg_data.groupby('genus')['recurrence'].sum()

    # Step 9: Create the stacked barplot with the custom color palette
    fig = px.bar(
        agg_data, y='genus', x='recurrence_normalized',
        title='Normalized Stacked Barplot of Predicted Superclasses Occurrence for Genus',
        labels={'recurrence_normalized': 'Normalized Recurrence (%)'},
        color='chemical_superclass',
        color_discrete_map=color_map,  # Use custom colors
        category_orders={'genus': unique_genera, 'chemical_superclass': unique_superclasses},
        orientation='h'
    )

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Genus<i>')

    # Set genus labels in italics
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=list(range(len(unique_genera))),
        ticktext=[f'<i>{genus}</i>' for genus in unique_genera]
    ))

    # Step 10: Add annotations for total compounds per genus
    max_x = agg_data['recurrence_normalized'].max()

    for genus, total_recurrence in total_recurrence_per_genus.items():
        # Position the annotation at the end of the bar with a slight offset
        x_position = max_x + 2
        if x_position > 100:  # Ensure the label is not outside the plot area
            x_position = 100  # Set it to 100 if it exceeds the maximum x-axis value

        # Add annotation to the plot
        fig.add_annotation(
            x=x_position,
            y=genus,
            text=f'Total compounds: {total_recurrence}',
            showarrow=False,
            font=dict(size=10, color='black'),
            xanchor='left',
            yanchor='middle'
        )

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Modify the figure size
    fig.update_layout(width=1500, height=1500)

    # Save the figure as an HTML file
    output_html_path = os.path.join(output_folder, 'Wikidata_superclass_barplot_genus_normalized.html')
    fig.write_html(output_html_path)

    # Show the figure
    fig.show()

    print(f"✅ Process completed! Visualization saved to {output_html_path}")

def plot_genus_pathway(output_folder):
    """
    Reads .tsv files, processes genus and pathway data, and generates a stacked bar plot.

    Parameters:
    - output_folder (str): Path to the folder containing 'genus_data' subfolder.
    - pathway_colors: Dictionary mapping pathways to specific colors.

    Saves the plot as an HTML file in the output folder.
    """
    # Define custom colors for the 7 pathway categories
    pathway_colors = {
        'Terpenoids': '#618264',  
        'Alkaloids': '#305F72',  
        'Shikimates and Phenylpropanoids': '#80558C',  
        'Polyketides': '#EF4B4B',  
        'Fatty acids': '#FF6C22',  
        'Amino acids and Peptides': '#F4E869',  
        'Carbohydrates': '#65451F'  
    } 

    genus_data_folder = os.path.join(output_folder, 'genus_data')

    # Step 1: Read data from all .tsv files
    all_data = pd.concat([
        pd.read_csv(os.path.join(genus_data_folder, filename), sep='\t')
        for filename in os.listdir(genus_data_folder) if filename.endswith(".tsv")
    ], ignore_index=True)

    # Step 2: Rename columns
    all_data.rename(columns={'organism_taxonomy_08genus': 'genus'}, inplace=True)
    all_data.rename(columns={'structure_taxonomy_npclassifier_01pathway': 'Pathway'}, inplace=True)

    # Step 3: Remove 'API Error' and 'Not Classified' from 'Pathway'
    all_data = all_data[~all_data['Pathway'].isin(['API Error', 'Not Classified'])]

    # Step 4: Group and aggregate data to calculate recurrence
    agg_data = all_data.groupby(['genus', 'Pathway']).size().reset_index(name='recurrence')

    # Sort the DataFrame by 'genus' alphabetically
    agg_data = agg_data.sort_values(by='genus')

    # Get unique genus names and pathways
    unique_genera = sorted(agg_data['genus'].unique())
    unique_pathways = sorted(agg_data['Pathway'].unique())

    # Step 5: Create the stacked barplot with custom colors
    fig = px.bar(
        agg_data, y=agg_data['genus'].apply(lambda x: f"<i>{x}</i>"), x='recurrence',
        title='Stacked Barplot of Predicted Pathways Occurrence for Genus',
        labels={'recurrence': 'Recurrence'},
        color='Pathway',
        color_discrete_map=pathway_colors,  # Use custom colors
        category_orders={'genus': unique_genera, 'Pathway': unique_pathways},  # Sort alphabetically
        orientation='h'
    )

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Genus<i>')

    # Set genus labels in italics
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=list(range(len(unique_genera))),
        ticktext=[f'<i>{genus}</i>' for genus in unique_genera]
    ))

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Modify the figure size
    fig.update_layout(width=1200, height=1400)

    # Save the figure as an HTML file
    output_html_path = os.path.join(output_folder, 'Wikidata_pathway_barplot_genus.html')
    fig.write_html(output_html_path)

    # Show the figure
    fig.show()

    print(f"✅ Process completed! Visualization saved to {output_html_path}")


def plot_genus_pathway_norm(output_folder):
    """
    Reads .tsv files, processes genus and pathway data, normalizes recurrence values, and generates a stacked bar plot.

    Parameters:
    - output_folder (str): Path to the folder containing 'genus_data' subfolder.
    - pathway_colors: Dictionary mapping pathways to specific colors.

    Saves the plot as an HTML file in the output folder.
    """
    # Define custom colors for the 7 pathway categories
    pathway_colors = {
        'Terpenoids': '#618264',  
        'Alkaloids': '#305F72',  
        'Shikimates and Phenylpropanoids': '#80558C',  
        'Polyketides': '#EF4B4B',  
        'Fatty acids': '#FF6C22',  
        'Amino acids and Peptides': '#F4E869',  
        'Carbohydrates': '#65451F'  
    } 
    genus_data_folder = os.path.join(output_folder, 'genus_data')

    # Step 1: Read data from all .tsv files
    all_data = pd.concat([
        pd.read_csv(os.path.join(genus_data_folder, filename), sep='\t')
        for filename in os.listdir(genus_data_folder) if filename.endswith(".tsv")
    ], ignore_index=True)

    # Step 2: Rename columns
    all_data.rename(columns={'organism_taxonomy_08genus': 'genus'}, inplace=True)
    all_data.rename(columns={'structure_taxonomy_npclassifier_01pathway': 'Pathway'}, inplace=True)

    # Step 3: Remove 'API Error' and 'Not Classified' from 'Pathway'
    all_data = all_data[~all_data['Pathway'].isin(['API Error', 'Not Classified'])]

    # Step 4: Group and aggregate data to calculate recurrence
    agg_data = all_data.groupby(['genus', 'Pathway']).size().reset_index(name='recurrence')

    # Step 5: Normalize the recurrence values within each genus group
    agg_data['recurrence_normalized'] = agg_data.groupby('genus')['recurrence'].transform(lambda x: x / x.sum()) * 100

    # Sort the DataFrame by 'genus' alphabetically
    agg_data = agg_data.sort_values(by='genus')

    # Get unique genus names and pathways
    unique_genera = sorted(agg_data['genus'].unique())
    unique_pathways = sorted(agg_data['Pathway'].unique())

    # Step 6: Calculate total recurrence for each genus
    total_recurrence_per_genus = agg_data.groupby('genus')['recurrence'].sum()

    # Step 7: Create the stacked barplot with custom colors
    fig = px.bar(
        agg_data, y=agg_data['genus'].apply(lambda x: f"<i>{x}</i>"), x='recurrence_normalized',
        title='Normalized Stacked Barplot of Predicted Pathways Occurrence for Genus',
        labels={'recurrence_normalized': 'Normalized Recurrence (%)'},
        color='Pathway',
        color_discrete_map=pathway_colors,  # Use custom colors
        category_orders={'genus': unique_genera, 'Pathway': unique_pathways},  # Sort alphabetically
        orientation='h'
    )

    # Modify the y-axis label
    fig.update_yaxes(title_text='<i>Genus<i>')

    # Set genus labels in italics
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=list(range(len(unique_genera))),
        ticktext=[f'<i>{genus}</i>' for genus in unique_genera]
    ))

    # Step 8: Add annotations for total compounds per genus
    max_x = agg_data['recurrence_normalized'].max()

    for genus, total_recurrence in total_recurrence_per_genus.items():
        # Position the annotation at the end of the bar with a slight offset
        x_position = max_x + 2
        if x_position > 100:  # Ensure the label is not outside the plot area
            x_position = 100  # Set it to 100 if it exceeds the maximum x-axis value

        # Add annotation to the plot
        fig.add_annotation(
            x=x_position,
            y=genus,
            text=f'Total compounds: {total_recurrence}',
            showarrow=False,
            font=dict(size=10, color='black'),
            xanchor='left',
            yanchor='middle'
        )

    # Set a white background
    fig.update_layout(plot_bgcolor='white')

    # Modify the figure size
    fig.update_layout(width=1200, height=1400)

    # Save the figure as an HTML file
    output_html_path = os.path.join(output_folder, 'Wikidata_pathway_barplot_genus_normalized.html')
    fig.write_html(output_html_path)

    # Show the figure
    fig.show()

    print(f"✅ Process completed! Visualization saved to {output_html_path}")


def dotplot_species_superclass(output_folder):
    """
    Reads .tsv files, processes species and superclass data, and generates a dot plot.

    Parameters:
    - output_folder (str): Path to the folder containing 'species_data' subfolder.

    Saves the plot as an HTML file in the output folder.
    """

    try:
        # Define base colors for each Pathway
        pathway_shades = {
            'Terpenoids': '#618264',
            'Alkaloids': '#305F72',
            'Shikimates and Phenylpropanoids': '#80558C',
            'Polyketides': '#EF4B4B',
            'Fatty acids': '#FF6C22',
            'Amino acids and Peptides': '#F4E869',
            'Carbohydrates': '#65451F'
        }

        species_data_folder = os.path.join(output_folder, 'species_data')

        # Step 1: Read data from all .tsv files
        all_data = pd.concat([
            pd.read_csv(os.path.join(species_data_folder, filename), sep='\t')
            for filename in os.listdir(species_data_folder) if filename.endswith(".tsv")
        ], ignore_index=True)

        # Step 2: Rename columns
        all_data.rename(columns={'organism_taxonomy_09species': 'species'}, inplace=True)

        # Function to split 'chemical_superclass' into pathway and superclass
        def split_chemical_superclass(row):
            if isinstance(row['chemical_superclass'], str):
                parts = row['chemical_superclass'].split('-')
                return parts if len(parts) == 2 else [parts[0], 'Unknown']
            return ['Unknown', 'Unknown']

        # Step 3: Apply the function to split chemical superclass
        all_data[['Pathway', 'superclass']] = all_data.apply(lambda row: split_chemical_superclass(row), axis=1, result_type='expand')

        # Step 4: Process data for color mapping
        color_map = {}
        for pathway, superclasses in all_data.groupby('Pathway')['chemical_superclass'].unique().items():
            for superclass in superclasses:
                color_map[f"{pathway}-{superclass}"] = pathway_shades.get(pathway, "#808080")  # Default to gray if missing

        # Step 5: Group and aggregate data to calculate recurrence
        agg_data = all_data.groupby(['species', 'Pathway', 'chemical_superclass']).size().reset_index(name='recurrence')

        # Convert 'species' column to categorical data
        agg_data['species'] = pd.Categorical(agg_data['species'], categories=sorted(agg_data['species'].unique()), ordered=True)

        # Get unique species names and sorted unique superclasses
        unique_species = sorted(agg_data['species'].unique())
        unique_superclasses = sorted(agg_data['chemical_superclass'].unique())

        # Step 6: Create a dot plot
        fig = px.scatter(
            agg_data, x='chemical_superclass', y='species', size='recurrence',
            labels={'chemical_superclass': 'Chemical Superclass', 'species': 'Species', 'recurrence': 'Recurrence'},
            title='Dot Plot of Recurrence of Chemical Superclasses for Species',
            color='Pathway',
            color_discrete_map=color_map,  # Use the color map
            category_orders={'chemical_superclass': unique_superclasses},
            width=1500, height=1500
        )

        # Modify the y-axis label
        fig.update_yaxes(title_text='<i>Species<i>')

        # Set species labels in italics
        fig.update_layout(yaxis=dict(
            tickmode='array',
            tickvals=list(range(len(unique_species))),
            ticktext=[f'<i>{species}</i>' for species in unique_species]
        ))

        # Set a white background
        fig.update_layout(plot_bgcolor='white')

        # Save the figure as an HTML file
        output_html_path = os.path.join(output_folder, 'Wikidata_sclass_dotplot_species.html')
        fig.write_html(output_html_path)

        # Show the figure
        fig.show()

        print(f"✅ Dot plot successfully created and saved: {output_html_path}")

    except Exception as e:
        print(f"❌ An error occurred: {str(e)}")

def dotplot_species_pathway(output_folder):
    """
    Reads .tsv files, processes species and pathway data, and generates a dot plot.

    Parameters:
    - output_folder (str): Path to the folder containing 'species_data' subfolder.

    Saves the plot as an HTML file in the output folder.
    """

    try:
        # Define fixed colors for each Pathway
        pathway_colors = {
            'Terpenoids': '#618264',
            'Alkaloids': '#305F72',
            'Shikimates and Phenylpropanoids': '#80558C',
            'Polyketides': '#EF4B4B',
            'Fatty acids': '#FF6C22',
            'Amino acids and Peptides': '#F4E869',
            'Carbohydrates': '#65451F'
        }

        species_data_folder = os.path.join(output_folder, 'species_data')

        # Step 1: Read data from all .tsv files
        all_data = pd.concat([
            pd.read_csv(os.path.join(species_data_folder, filename), sep='\t')
            for filename in os.listdir(species_data_folder) if filename.endswith(".tsv")
        ], ignore_index=True)

        # Step 2: Rename columns
        all_data.rename(columns={'organism_taxonomy_09species': 'species'}, inplace=True)
        all_data.rename(columns={'structure_taxonomy_npclassifier_01pathway': 'Pathway'}, inplace=True)

        # Step 3: Remove 'API Error' and 'Not Classified' from 'Pathway'
        all_data = all_data[~all_data['Pathway'].isin(['API Error', 'Not Classified'])]

        # Step 4: Group and aggregate data to calculate recurrence
        agg_data = all_data.groupby(['species', 'Pathway']).size().reset_index(name='recurrence')

        # Convert 'species' column to categorical data
        agg_data['species'] = pd.Categorical(agg_data['species'], categories=sorted(agg_data['species'].unique()), ordered=True)

        # Get unique species names and pathways
        unique_species = sorted(agg_data['species'].unique())
        unique_pathways = sorted(agg_data['Pathway'].unique())

        # Step 5: Create a dot plot
        fig = px.scatter(
            agg_data, x='Pathway', y='species', size='recurrence',
            labels={'Pathway': 'Pathway', 'species': 'Species', 'recurrence': 'Recurrence'},
            title='Dot Plot of Recurrence of Pathways for Species',
            color='Pathway',
            color_discrete_map=pathway_colors,
            width=1500, height=3000
        )

        # Modify the y-axis label
        fig.update_yaxes(title_text='<i>Species<i>')

        # Set species labels in italics
        fig.update_layout(yaxis=dict(
            tickmode='array',
            tickvals=list(range(len(unique_species))),
            ticktext=[f'<i>{species}</i>' for species in unique_species]
        ))

        # Set a white background
        fig.update_layout(plot_bgcolor='white')

        # Save the figure as an HTML file
        output_html_path = os.path.join(output_folder, 'Wikidata_pathway_dotplot_species.html')
        fig.write_html(output_html_path)

        # Show the figure
        fig.show()

        print(f"✅ Dot plot successfully created and saved: {output_html_path}")

    except Exception as e:
        print(f"❌ An error occurred: {str(e)}")


def heatmap_pathway_species(output_folder):
    """
    Reads .tsv files, processes species and pathway data, normalizes recurrence values, and generates a heatmap.

    Parameters:
    - output_folder (str): Path to the folder containing 'species_data' subfolder.

    Saves the heatmap as an HTML file in the output folder.
    """

    try:
        species_data_folder = os.path.join(output_folder, 'species_data')

        # Step 1: Read data from all .tsv files
        all_data = pd.concat([
            pd.read_csv(os.path.join(species_data_folder, filename), sep='\t')
            for filename in os.listdir(species_data_folder) if filename.endswith(".tsv")
        ], ignore_index=True)

        # Step 2: Rename columns
        all_data.rename(columns={'organism_taxonomy_09species': 'species'}, inplace=True)
        all_data.rename(columns={'structure_taxonomy_npclassifier_01pathway': 'Pathway'}, inplace=True)

        # Step 3: Remove 'API Error' and 'Not Classified' from 'Pathway'
        all_data = all_data[~all_data['Pathway'].isin(['API Error', 'Not Classified'])]

        # Step 4: Group and aggregate data to calculate recurrence
        agg_data = all_data.groupby(['species', 'Pathway']).size().reset_index(name='recurrence')

        # Step 5: Normalize the recurrence values within each species group
        agg_data['recurrence_normalized'] = agg_data.groupby('species')['recurrence'].transform(lambda x: x / x.sum()) * 100

        # Convert 'species' column to categorical data
        agg_data['species'] = pd.Categorical(agg_data['species'], categories=sorted(agg_data['species'].unique()), ordered=True)

        # Get unique species names and pathways
        unique_species = sorted(agg_data['species'].unique())
        unique_pathways = sorted(agg_data['Pathway'].unique())

        # Step 6: Create a DataFrame with all possible combinations of species and pathways
        all_combinations = pd.DataFrame([(species, pathway) for species in unique_species for pathway in unique_pathways],
                                        columns=['species', 'Pathway'])

        # Merge the all_combinations DataFrame with the aggregated data to include all species
        merged_data = pd.merge(all_combinations, agg_data, on=['species', 'Pathway'], how='left')

        # Fill missing values (NaN) with 0 for recurrence and recurrence_normalized columns
        merged_data['recurrence'].fillna(0, inplace=True)
        merged_data['recurrence_normalized'].fillna(0, inplace=True)

        # Step 7: Pivot the merged data to have 'species' as rows and 'Pathway' as columns
        pivot_data = merged_data.pivot_table(index='species', columns='Pathway', values='recurrence_normalized', fill_value=0)

        # Step 8: Create the heatmap
        fig = px.imshow(
            pivot_data,
            labels=dict(x="Pathway", y="Species", color="Normalized Recurrence (%)"),
            x=pivot_data.columns,
            y=pivot_data.index,
            color_continuous_scale='tealrose',  # You can choose any other color scale
            title='Heatmap of Normalized Recurrence of Pathways for Species'
        )

        # Modify the y-axis label
        fig.update_yaxes(title_text='<i>Species<i>')

        # Set species labels in italics
        fig.update_layout(yaxis=dict(
            tickmode='array',
            tickvals=list(range(len(unique_species))),
            ticktext=[f'<i>{species}</i>' for species in unique_species]
        ))

        # Modify the size of the figure
        fig.update_layout(width=1500, height=1500)

        # Save the figure as an HTML file
        output_html_path = os.path.join(output_folder, 'Wikidata_pathway_heatmap_species_normalized.html')
        fig.write_html(output_html_path)

        # Show the figure
        fig.show()

        print(f"✅ Heatmap successfully created and saved: {output_html_path}")

    except Exception as e:
        print(f"❌ An error occurred: {str(e)}")

