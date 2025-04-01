from rcbgk import *

cpdag_input = '{"0": {"1", "2"}, "1": {"0", "2"}, "2": {"0", "1"}}'
bgk_input = '[("0", "1", 1)]'
cpdag_input = '{"A": ["J", "C", "D", "E", "I"], "B": ["F", "G"], "C": ["A", "I"], "D": ["A", "E", "H", "F"], "E": ["D", "A", "I"], "F": ["G"], "G": [], "H": ["D"], "I": [], "J": ["A"]}'
bgk_input = '[("J","D",1)]'

input_dict = eval(cpdag_input)
bgk = eval(bgk_input) if bgk_input else []

# Create a mapping between original labels and numeric indices
nodes = sorted(input_dict.keys())
label_to_index = {node: i for i, node in enumerate(nodes)}
index_to_label = {i: node for node, i in label_to_index.items()}

# Convert dictionary to NetworkX DiGraph
input_graph = nx.DiGraph()
input_graph.add_nodes_from(range(len(nodes)))

# Add edges maintaining CPDAG structure
for node, neighbors in input_dict.items():
    node_idx = label_to_index[node]
    for neighbor in neighbors:
        neighbor_idx = label_to_index[neighbor]
        # Only add edge in one direction for directed edges
        if node in input_dict.get(neighbor, []):
            # Undirected edge - add both directions
            input_graph.add_edge(node_idx, neighbor_idx)
            input_graph.add_edge(neighbor_idx, node_idx)
        else:
            # Directed edge - add only one direction
            input_graph.add_edge(node_idx, neighbor_idx)

# Convert background knowledge to numeric indices
numeric_bgk = []
for x, y, t in bgk:
    x_idx = label_to_index[x] if isinstance(x, str) else x
    y_idx = label_to_index[y] if isinstance(y, str) else y
    numeric_bgk.append((x_idx, y_idx, t))

# Convert CPDAG and background knowledge to MPDAG
dcc = transCbgk(input_graph, numeric_bgk)
mpdag = construct_mpdag(input_graph, dcc)
print(f"dcc: {dcc}")
print(f"mpdag: {mpdag}")

# Convert back to original labels
output_graph = nx.DiGraph()
output_graph.add_nodes_from(nodes)

# Add edges while preserving direction
for u, v in mpdag.edges():
    u_label = index_to_label[u]
    v_label = index_to_label[v]
    output_graph.add_edge(u_label, v_label)

print(output_graph.edges)

# Create visualization
image_data = plot_graph(output_graph, title="MPDAG with Background Knowledge", output_path="./test_mpdag.png")