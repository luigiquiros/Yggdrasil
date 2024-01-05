import opentree
import ete4

def main():
    try:
        tree = opentree.OT.taxon_subtree(ott_id='601168')  # ott_id for Celastraceae
    except Exception as e:
        print(f"Failed to access Open Tree of Life API: {e}")
        return

    try:
        newick_str = tree.response_dict['newick']
        t = ete4.Tree(newick_str)
        t.show()
    except Exception as e:
        print(f"Error parsing Newick string: {e}")
        return

if __name__ == "__main__":
    main()


from ete4 import Tree

t = Tree(open('Celastraceae.tre'))

t.explore()

### celastraceae tree downloaded and smartview

from ete4 import Tree
from ete4.smartview import TreeLayout, RectFace, TextFace, ScaleFace
import random

t = Tree(open('Celastraceae.tre'))

# annotate numerical values to each leaf
for node in t.leaves():
    node.add_prop('count', random.randint(1, 100))

# define tree style function
def layout_tree_style(tree_style):
    # add scale bar to footer
    scaleface = ScaleFace(
        name='sample1', 
        width=150, 
        color='black',
        scale_range=(0, 100), 
        tick_width=80, 
        line_width=1,
        formatter='%.0f',
        min_fsize=6, 
        max_fsize=12, 
        ftype='sans-serif',
        padding_x=0, 
        padding_y=0)

    tree_style.aligned_panel_header.add_face(scaleface, column=0)
    tree_style.aligned_panel_footer.add_face(scaleface, column=0)

    # add title to header and footer
    text = TextFace("Count", min_fsize=5, max_fsize=12, width=50, rotation=0)
    tree_style.aligned_panel_header.add_face(text, column=0)    
    return 

# define node Face layout function
def layout_barplot(node):
    if node.is_leaf:
        width = node.props.get('count') * 1.5
        rect_face = RectFace(
            width=width, height=70, color='skyblue',
            opacity=0.7, text=None, fgcolor='black',
            min_fsize=6, max_fsize=15, ftype='sans-serif',
            padding_x=0, padding_y=0,
            tooltip=None)
        node.add_face(rect_face, position='aligned', column=0)
        return 

# Create a TreeLayout object, passing in the function
barplot_layout = TreeLayout(
    name='BarPlot',
    ns=layout_barplot, 
    ts=layout_tree_style,
    aligned_faces=True)

# add layout to layouts list
layouts = []
layouts.append(barplot_layout)
t.explore(
    layouts=layouts, 
    include_props=("name", "dist", "count"),
    keep_server=True)



##### test recovery from a list of sp ######

import csv
import requests
import time

input_path = '/mnt/c/Users/quirosgu/Desktop/Inmuno/Clean_collection_taxonomical_data.csv'
species_header = 'query_otol_species'

# Configure the rate limit
REQUESTS_PER_MINUTE = 10
REQUEST_INTERVAL = 60.0 / REQUESTS_PER_MINUTE  # time between requests in seconds

def get_ott_id(species_name):
    time.sleep(REQUEST_INTERVAL)  # Wait before making the request
    response = requests.get(f"https://api.opentreeoflife.org/v3/tnrs/match_names", json={'names': [species_name]})
    if response.status_code == 200:
        results = response.json()
        return results['results'][0]['matches'][0]['taxon']['ott_id'] if results['results'][0]['matches'] else None
    else:
        print(f"Error fetching OTT ID for {species_name}")
        return None

def get_newick_tree(ott_ids):
    time.sleep(REQUEST_INTERVAL)  # Wait before making the request
    response = requests.post("https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree", json={'ott_ids': ott_ids})
    if response.status_code == 200:
        return response.json()['newick']
    else:
        print("Error fetching Newick tree")
        return None

def read_species_from_csv(file_path):
    species_list = []
    with open(file_path, mode='r', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            species_list.append(row[species_header])
    return species_list

def main():
    file_path = input_path # Replace with your CSV file path
    species_list = read_species_from_csv(file_path)

    ott_ids = []
    for species in species_list:
        ott_id = get_ott_id(species)
        if ott_id:
            ott_ids.append(ott_id)

    if ott_ids:
        newick_tree = get_newick_tree(ott_ids)
        if newick_tree:
            t = Tree(newick_tree)
            # annotate numerical values to each leaf
            for node in t.leaves():
                node.add_prop('count', random.randint(1, 100))
            # define tree style function
            
            def layout_tree_style(tree_style):
                # add scale bar to footer
                scaleface = ScaleFace(
                    name='sample1', 
                    width=150, 
                    color='black',
                    scale_range=(0, 100), 
                    tick_width=80, 
                    line_width=1,
                    formatter='%.0f',
                    min_fsize=6, 
                    max_fsize=12, 
                    ftype='sans-serif',
                    padding_x=0, 
                    padding_y=0)

                tree_style.aligned_panel_header.add_face(scaleface, column=0)
                tree_style.aligned_panel_footer.add_face(scaleface, column=0)

                # add title to header and footer
                text = TextFace("Count", min_fsize=5, max_fsize=12, width=50, rotation=0)
                tree_style.aligned_panel_header.add_face(text, column=0)    
                return 
            
            # define node Face layout function
            def layout_barplot(node):
                if node.is_leaf:
                    width = node.props.get('count') * 1.5
                    rect_face = RectFace(
                        width=width, height=70, color='skyblue',
                        opacity=0.7, text=None, fgcolor='black',
                        min_fsize=6, max_fsize=15, ftype='sans-serif',
                        padding_x=0, padding_y=0,
                        tooltip=None)
                    node.add_face(rect_face, position='aligned', column=0)
                    return 

            # Create a TreeLayout object, passing in the function
            barplot_layout = TreeLayout(
                name='BarPlot',
                ns=layout_barplot, 
                ts=layout_tree_style,
                aligned_faces=True)

            # add layout to layouts list
            layouts = []
            layouts.append(barplot_layout)
            t.explore(
                layouts=layouts, 
                include_props=("name", "dist", "count"),
                keep_server=True)
                        
        else:
            print("Failed to retrieve Newick tree.")
    else:
        print("Failed to resolve OTT IDs for all species.")

if __name__ == "__main__":
    main()


### opt 2

import csv
import requests
import time
import random
from ete4 import Tree
from ete4.smartview import TreeLayout, RectFace, TextFace, ScaleFace, TreeStyle
from tqdm import tqdm 

input_path = '/mnt/c/Users/quirosgu/Desktop/Celastraceae/Celastraceae1.csv'#Clean_collection_taxonomical_data.csv'
species_header = 'ATTRIBUTE_Species'#'''query_otol_species'

REQUESTS_PER_MINUTE = 30
REQUEST_INTERVAL = 60.0 / REQUESTS_PER_MINUTE

def get_ott_id(species_name):
    """
    Retrieves the OTT ID for a given species name using the Open Tree of Life API.
    Args:
        species_name (str): The name of the species.
    Returns:
        int: The OTT ID of the species, or None if not found or in case of an error.
    """
    time.sleep(REQUEST_INTERVAL)
    try:
        response = requests.post("https://api.opentreeoflife.org/v3/tnrs/match_names", json={'names': [species_name]})
        if response.status_code == 200:
            results = response.json()
            print(f"Response for {species_name}: {results}")  # Debugging: Print the response
            return results['results'][0]['matches'][0]['taxon']['ott_id'] if results['results'][0]['matches'] else None
        else:
            print(f"Error fetching OTT ID for {species_name}: HTTP {response.status_code}, Response: {response.text}")
            return None
    except requests.RequestException as e:
        print(f"Request exception occurred for {species_name}: {e}")
        return None


def get_newick_tree(ott_ids):
    """
    Retrieves the Newick tree for a list of OTT IDs using the Open Tree of Life API.
    Args:
        ott_ids (list of int): The OTT IDs for which to retrieve the tree.
    Returns:
        str: The Newick tree string, or None in case of an error.
    """
    time.sleep(REQUEST_INTERVAL)
    try:
        response = requests.post("https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree", json={'ott_ids': ott_ids})
        if response.status_code == 200:
            return response.json()['newick']
        else:
            print("Error fetching Newick tree: HTTP " + str(response.status_code))
            return None
    except requests.RequestException as e:
        print(f"Request exception occurred: {e}")
        return None

def read_species_from_csv(file_path):
    """
    Reads species names from a CSV file.
    Args:
        file_path (str): The path to the CSV file.
    Returns:
        list of str: A list of species names.
    """
    species_list = []
    with open(file_path, mode='r', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            species_list.append(row[species_header])
    return species_list

def main():
    species_list = read_species_from_csv(input_path)

    # Initialize the progress bar
    pbar = tqdm(total=len(species_list), desc="Processing species", unit="species")

    ott_ids = []
    for species in species_list:
        ott_id = get_ott_id(species)
        if ott_id:
            ott_ids.append(ott_id)
        pbar.update(1)  # Update the progress bar for each species processed

    pbar.close()  # Close the progress bar

    if ott_ids:
        newick_tree = get_newick_tree(ott_ids)
        if newick_tree:
            t = Tree(newick_tree)
            
            # Create a Tree Style
            ts = TreeStyle()
            ts.show_leaf_name = True
            ts.show_branch_length = True
            ts.show_branch_support = True
            ts.scale = 120 # You can adjust this scale factor
            ts.mode = "c"  # Use circular tree layout

            # Render the tree and save it to a file
            t.render("tree_visualization.png", w=500, tree_style=ts)  # Adjust the width as needed

            # Render the tree to HTML
            tree_html = t.write()  # Get Newick string with internal node names
            html_output = f"<html><body><pre>{tree_html}</pre></body></html>"

            # Save the HTML to a file
            with open("tree_visualization.html", "w") as file:
                file.write(html_output)

            print("Tree visualization saved as HTML.")
        else:
            print("Failed to retrieve Newick tree.")
    else:
        print("Failed to resolve OTT IDs for all species.")
        
if __name__ == "__main__":
    main()
