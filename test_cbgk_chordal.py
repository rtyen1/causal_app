import networkx as nx
import numpy as np
from rcbgk import (genConsistentCbgk, transCbgk, is_consistent, construct_mpdag, 
                  generate_undirected_chordal, dag_to_cpdag, bgk_ida, true_cov)

def test_with_bgk():
    # Create a 10-node chordal graph with edge ratio 0.3
    n_nodes = 20
    edge_ratio = 0.4
    
    # Generate undirected chordal graph
    chordal_graph = generate_undirected_chordal(n_nodes, edge_ratio)
    print("Generated chordal graph edges:", list(chordal_graph.edges()))
    
    # Convert to DAG by giving random orientations that maintain acyclicity
    dag = nx.DiGraph()
    dag.add_nodes_from(chordal_graph.nodes())
    
    # Add edges with orientations that follow alphabetical order to ensure acyclicity
    for u, v in chordal_graph.edges():
        if u < v:
            dag.add_edge(u, v)
            dag[u][v]['weight'] = np.random.uniform(0.3, 0.7)
        else:
            dag.add_edge(v, u)
            dag[v][u]['weight'] = np.random.uniform(0.3, 0.7)
            
    print("\nDAG edges with weights:", [(u, v, d['weight']) for u, v, d in dag.edges(data=True)])
    
    # Convert to CPDAG
    cpdag = dag_to_cpdag(dag)
    print("\nCPDAG edges:", list(cpdag.edges()))
    
    # Calculate covariance matrix
    mcov = true_cov(dag)
    print("\nCovariance matrix:")
    print(mcov)
    
    # Sample x and y nodes
    nodes = list(dag.nodes())
    x, y = np.random.choice(nodes, size=2, replace=False)
    print(f"\nSelected nodes for causal effect estimation: x={x}, y={y}")
    
    # Generate three types of background knowledge
    num_constraints = 2  # Generate 2 constraints of each type
    
    # Direct causal constraints
    bgk_d = genConsistentCbgk(dag, num_constraints, bgk_type=1)
    print("\nDirect causal background knowledge:", bgk_d)
    dcc_d = transCbgk(cpdag, bgk_d)
    print("Direct causal DCC:", dcc_d)
    
    # Non-causal constraints
    bgk_nc = genConsistentCbgk(dag, num_constraints, bgk_type=2)
    print("\nNon-causal background knowledge:", bgk_nc)
    dcc_nc = transCbgk(cpdag, bgk_nc)
    print("Non-causal DCC:", dcc_nc)
    
    # Causal constraints
    bgk_c = genConsistentCbgk(dag, num_constraints, bgk_type=3)
    print("\nCausal background knowledge:", bgk_c)
    dcc_c = transCbgk(cpdag, bgk_c)
    print("Causal DCC:", dcc_c)
    
    # Check consistency for each type
    print("\nConsistency check:")
    
    is_consistent_d = is_consistent(cpdag, dcc_d)
    print("Direct causal consistency:", is_consistent_d['indicator'])
    if is_consistent_d['indicator']:
        print("Direct causal PEO:", is_consistent_d['PEO'])
        
    is_consistent_nc = is_consistent(cpdag, dcc_nc)
    print("\nNon-causal consistency:", is_consistent_nc['indicator'])
    if is_consistent_nc['indicator']:
        print("Non-causal PEO:", is_consistent_nc['PEO'])
        
    is_consistent_c = is_consistent(cpdag, dcc_c)
    print("\nCausal consistency:", is_consistent_c['indicator'])
    if is_consistent_c['indicator']:
        print("Causal PEO:", is_consistent_c['PEO'])
    
    # Calculate causal effects with different types of background knowledge
    print("\nCalculating causal effects:")
    
    # Without background knowledge
    effects_no_bgk = bgk_ida(x, y, mcov, dag)
    print(f"\nPossible effects without background knowledge: {effects_no_bgk}")
    
    # With direct causal background knowledge
    if is_consistent_d['indicator']:
        effects_d = bgk_ida(x, y, mcov, dag, bgk_d)
        print(f"Possible effects with direct causal background knowledge: {effects_d}")
    
    # With non-causal background knowledge
    if is_consistent_nc['indicator']:
        effects_nc = bgk_ida(x, y, mcov, dag, bgk_nc)
        print(f"Possible effects with non-causal background knowledge: {effects_nc}")
    
    # With causal background knowledge
    if is_consistent_c['indicator']:
        effects_c = bgk_ida(x, y, mcov, dag, bgk_c)
        print(f"Possible effects with causal background knowledge: {effects_c}")
    
    # If all consistent, construct MPDAG for each type
    if all([is_consistent_d['indicator'], is_consistent_nc['indicator'], is_consistent_c['indicator']]):
        print("\nConstructing MPDAGs:")
        
        mpdag_d = construct_mpdag(cpdag, dcc_d)
        print("\nDirect causal MPDAG edges:", list(mpdag_d.edges()))
        
        mpdag_nc = construct_mpdag(cpdag, dcc_nc)
        print("\nNon-causal MPDAG edges:", list(mpdag_nc.edges()))
        
        mpdag_c = construct_mpdag(cpdag, dcc_c)
        print("\nCausal MPDAG edges:", list(mpdag_c.edges()))
        
        # Print statistics
        print("\nStatistics:")
        print(f"Number of edges in original DAG: {len(dag.edges())}")
        print(f"Number of edges in CPDAG: {len(cpdag.edges())}")
        print(f"Number of edges in direct causal MPDAG: {len(mpdag_d.edges())}")
        print(f"Number of edges in non-causal MPDAG: {len(mpdag_nc.edges())}")
        print(f"Number of edges in causal MPDAG: {len(mpdag_c.edges())}")
        
        # Print effect size statistics
        print("\nEffect size statistics:")
        print(f"Number of possible effects without background knowledge: {len(effects_no_bgk)}")
        if is_consistent_d['indicator']:
            print(f"Number of possible effects with direct causal background knowledge: {len(effects_d)}")
        if is_consistent_nc['indicator']:
            print(f"Number of possible effects with non-causal background knowledge: {len(effects_nc)}")
        if is_consistent_c['indicator']:
            print(f"Number of possible effects with causal background knowledge: {len(effects_c)}")

if __name__ == "__main__":
    test_with_bgk() 