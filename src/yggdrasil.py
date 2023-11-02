import opentreelib
import plotly.graph_objs as go

# Access Open Tree of Life API
tree = opentreelib.get_tree_of_life()
# Replace 'Celastraceae' with the desired taxonomic name or ID
taxon = tree.taxon("Celastraceae")
subtree = taxon.subtree()

# Process the Taxonomic Tree
# Example: Extract and prepare data for Plotly
# Assuming 'subtree' contains the necessary taxonomic data
# Process the tree structure to extract node names, parent-child relationships, etc.

node_labels = []  # Contains node names for the plot
parent_indices = []  # Contains indices of parent nodes
child_indices = []  # Contains indices of child nodes

# Use the 'subtree' data to populate node_labels, parent_indices, child_indices

# Visualize with Plotly
fig = go.Figure(go.Sankey(
    node=dict(
        label=node_labels,
    ),
    link=dict(
        source=parent_indices,
        target=child_indices
    )
))

fig.update_layout(title_text="Taxonomic Tree of Celastraceae Family")
fig.show()
