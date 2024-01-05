import dash
import dash_cytoscape as cyto
from dash import html
import requests
from Bio import Phylo
from io import StringIO

def fetch_and_parse_tree(family_name):
    api_url = "https://api.opentreeoflife.org/v3/tree_of_life/subtree"
    payload = {"ott_id": family_name, "format": "newick"}

    try:
        response = requests.post(api_url, json=payload)
        response.raise_for_status()  
        tree_data = response.json()["newick"]
        tree = Phylo.read(StringIO(tree_data), "newick")
        return tree
    except requests.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except Exception as err:
        print(f"An error occurred: {err}")

def convert_tree_for_dash_phylogeny(tree):
    # Convert the tree to a format suitable for dash-phylogeny
    pass

app = dash.Dash(__name__)

family_name = "Hominidae"  # Example family name
tree = fetch_and_parse_tree(family_name)
tree_data_for_dash_phylogeny = convert_tree_for_dash_phylogeny(tree)

app.layout = html.Div([
    dash_phylogeny.Tree(  # Assuming Tree is the component name
        id='phylogeny-tree',
        data=tree_data_for_dash_phylogeny,
        # Include other necessary properties as required by dash-phylogeny
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)