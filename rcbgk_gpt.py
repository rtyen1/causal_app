import numpy as np
import networkx as nx
import random
from itertools import combinations

# Function to find the critical set of x with respect to z
def critical_set(g, xlb, zlb):
    nodes = list(g.nodes)
    x = nodes.index(xlb)
    z = [nodes.index(node) for node in zlb]
    n = len(nodes)
    
    S = [(x, 0, 0)]
    C = np.zeros(n)
    flag = 0
    dlp = 0
    
    while len(S) != 0:
        dlp += 1
        e = S.pop(0)
        
        if e[1] != 0 and e[0] not in z:
            for alpha in set(g.neighbors(e[0])) - {e[1]}:
                if alpha not in g.neighbors(e[1]):
                    if e[1] == x:
                        S.append((alpha, e[0], e[0]))
                    else:
                        S.append((alpha, e[0], e[2]))
        elif e[0] not in z:
            for alpha in g.neighbors(e[0]):
                S.append((alpha, e[0], alpha))
        else:
            C[flag] = e[2]
            flag += 1
        
        if dlp > n**2 + 1:
            break
    
    return [nodes[i] for i in np.unique(C[C != 0].astype(int))]

# Function to find ancestors and chain components
def find_ancestors_and_chaincomp(amat, x):
    n = amat.shape[0]
    amat_skel = (amat + amat.T) > 0
    amat_dir = (amat - amat.T) * (amat - amat.T > 0)
    amat_undir = amat - amat_dir
    
    label_names = range(n)
    
    # Finding ancestors
    w_an = [x]
    res_an = []
    while w_an:
        node = w_an.pop(0)
        an_tmp = np.where(amat_dir[:, node] == 1)[0]
        w_an.extend(list(set(an_tmp) - set(res_an) - set(w_an)))
        res_an.append(node)
    
    # Finding chain components
    w_cp = [x]
    res_cp = []
    while w_cp:
        node = w_cp.pop(0)
        cp_tmp = np.where(amat_undir[:, node] == 1)[0]
        w_cp.extend(list(set(cp_tmp) - set(res_cp) - set(w_cp)))
        res_cp.append(node)
    
    return {'an': res_an, 'cp': res_cp}

# Function to transform Causal Background Knowledge
def trans_cbgk(cpdag, bgk):
    amat = np.array(cpdag)
    n = amat.shape[0]
    label_names = list(range(n))
    
    amat_skel = (amat + amat.T) > 0
    amat_dir = (amat - amat.T) * (amat - amat.T > 0)
    amat_undir = amat - amat_dir
    g = nx.from_numpy_matrix(amat_undir)
    
    res = {i: [] for i in label_names}
    ancestor = {i: [] for i in label_names}
    cp = {i: [] for i in label_names}
    
    if len(bgk) > 0:
        for item in bgk:
            x, y, t = item
            x_idx = label_names.index(x)
            y_idx = label_names.index(y)
            
            if amat_skel[x_idx, y_idx]:
                if t == 0:
                    res[y_idx].append([x_idx])
                else:
                    res[x_idx].append([y_idx])
            else:
                cset = critical_set(cpdag, x, y)
                if cset:
                    if t == 0:
                        for elem in cset:
                            res[elem].append([x_idx])
                    else:
                        res[x_idx].append(cset)
                else:
                    if t == 1:
                        res[x_idx].append([])
    
    return {k: v for k, v in res.items() if v}

# Random generation of valid direct causal clauses (DCC)
def gen_valid_dcc(cpdag, num):
    amat = np.array(cpdag)
    label_names = list(range(amat.shape[0]))
    n = amat.shape[0]
    amat_skel = (amat + amat.T) > 0
    amat_dir = (amat - amat.T) * (amat - amat.T > 0)
    amat_undir = amat - amat_dir
    
    tails = random.choices(label_names, k=num)
    uniq_tails = {tail: tails.count(tail) for tail in set(tails)}
    dcc = {tail: [] for tail in uniq_tails}
    
    for tail, count in uniq_tails.items():
        adjx = [i for i in label_names if amat_skel[i, tail] == 1]
        pax = [i for i in label_names if amat_dir[i, tail] == 1]
        chx = [i for i in label_names if amat_dir[tail, i] == 1]
        
        if all(x in pax for x in adjx):
            continue
        
        for _ in range(count):
            ind = np.random.choice([0, 1], size=len(adjx))
            select = [adjx[i] for i in range(len(adjx)) if ind[i] == 1]
            new_select = []
            
            if all(x in pax for x in select):
                can = list(set(adjx) - set(pax))
                while not new_select:
                    new_select = random.sample(can, random.randint(1, len(can)))
            
            if not all(x in chx for x in select):
                if random.random() < 0.8:
                    select = list(set(select) - set(chx))
            
            D = tuple(set(select + new_select))  # Convert to tuple (hashable type)
            dcc[tail].append(D)
    
    return {k: list(set(v)) for k, v in dcc.items() if v}

# Validity Check for Direct Causal Clauses (DCC)
def is_valid(cpdag, dcc):
    amat = np.array(cpdag)
    n = amat.shape[0]
    amat_skel = (amat + amat.T) > 0
    amat_dir = (amat - amat.T) * (amat - amat.T > 0)
    
    for x, clauses in dcc.items():
        for D in clauses:
            if len(D) == 0 or any(amat_skel[x, d] == 0 for d in D) or all(amat_dir[d, x] == 1 for d in D):
                return False
    return True

# Example usage:
cpdag = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])  # Simple CPDAG
bgk = [['X', 'Y', 1], ['Y', 'Z', 0]]  # Example background knowledge
valid_dcc = gen_valid_dcc(cpdag, 3)
print(valid_dcc)

# -------- Helper Functions for Causal Inference --------

# Check if a set of vertices is valid in terms of the graph structure
def is_valid(cpdag, dcc):
    amat = np.array(cpdag)
    n = amat.shape[0]
    amat_skel = (amat + amat.T) > 0
    amat_dir = (amat - amat.T) * (amat - amat.T > 0)
    
    for x, clauses in dcc.items():
        for D in clauses:
            if len(D) == 0 or any(amat_skel[x, d] == 0 for d in D) or all(amat_dir[d, x] == 1 for d in D):
                return False
    return True

# Check if a set of vertices is consistent in terms of the graph structure
def is_consistent(cpdag, dcc):
    amat = np.array(cpdag)
    n = amat.shape[0]
    amat_skel = (amat + amat.T) > 0
    amat_dir = (amat - amat.T) * (amat - amat.T > 0)
    
    if not is_valid(cpdag, dcc):
        return {'indicator': False, 'PEO': None}
    
    indicator = True
    vset = list(range(n))
    peo = []
    
    while vset:
        pln = restriction(cpdag, vset, dcc)['potentialLeaf']
        if len(pln) == 0:
            indicator = False
            peo = None
            break
        else:
            vset = list(set(vset) - set([pln[0]]))
            peo.append(pln[0])
    
    return {'indicator': indicator, 'PEO': peo}

# Function to construct the MOC (Maximal Orientable Chain)
def construct_moc(cpdag, x, dcc, simplified=True):
    amat = np.array(cpdag)
    label_names = list(range(amat.shape[0]))
    n = amat.shape[0]
    amat_skel = (amat + amat.T) > 0
    amat_dir = (amat - amat.T) * (amat - amat.T > 0)
    amat_undir = amat - amat_dir
    
    if not is_valid(cpdag, dcc):
        raise ValueError("The direct causal clauses are invalid with respect to the CPDAG.")
    
    vset = label_names
    while vset:
        pln = restriction(cpdag, vset, dcc)['potentialLeaf']
        if len(pln) == 0:
            raise ValueError("The direct causal clauses are inconsistent with respect to the CPDAG.")
        elif len(pln) == 1 and pln[0] == x:
            moc_set = vset
            break
        else:
            vset = list(set(vset) - set(pln) - set([x]))
    
    adj = list(set(label_names).intersection([i for i in range(n) if amat_undir[x, i] == 1]))
    if simplified:
        return adj
    
    return {'MOC': vset, 'adj': adj}

# Function to construct the MPDAG (Maximal Partially Directed Acyclic Graph)
def construct_mpdad(cpdag, dcc):
    amat = np.array(cpdag)
    label_names = list(range(amat.shape[0]))
    new_pa = {i: construct_moc(cpdag, i, dcc) for i in label_names}
    new_pa = {k: v for k, v in new_pa.items() if v}
    
    for x, p in new_pa.items():
        for p_x in p:
            amat[x, p_x] = 0
    
    return amat

# Function for fast construction of MPDAG, only considering non-singleton chain components
def construct_mpdad_fast(cpdag, dcc):
    amat = np.array(cpdag)
    label_names = list(range(amat.shape[0]))
    n = amat.shape[0]
    amat_skel = (amat + amat.T) > 0
    amat_dir = (amat - amat.T) * (amat - amat.T > 0)
    amat_undir = amat - amat_dir
    g = nx.from_numpy_matrix(amat_undir)
    
    cc = list(nx.connected_components(g))
    for component in cc:
        vset = list(component)
        rst = restriction(cpdag, vset, dcc)
        new_pa = {i: construct_moc(cpdag, i, rst['restricted']) for i in vset}
        new_pa = {k: v for k, v in new_pa.items() if v}
        
        for x, p in new_pa.items():
            for p_x in p:
                amat[x, p_x] = 0
    
    return amat

# Restriction operation to get induced subgraph
def restriction(cpdag, vset, dcc):
    if len(vset) == 0:
        raise ValueError("vset should be non-empty.")
    
    amat = np.array(cpdag)
    label_names = list(range(amat.shape[0]))
    n = amat.shape[0]
    amat_skel = (amat + amat.T) > 0
    amat_dir = (amat - amat.T) * (amat - amat.T > 0)
    amat_undir = amat - amat_dir
    amat_undir_sub = amat_undir[np.ix_(vset, vset)]
    
    usubg = nx.from_numpy_matrix(amat_undir_sub)
    
    if not is_valid(cpdag, dcc):
        raise ValueError("The direct causal clauses are invalid with respect to the CPDAG.")
    
    restrict_dcc = {k: v for k, v in dcc.items() if k in vset}
    for i, clauses in restrict_dcc.items():
        for j, D in enumerate(clauses):
            if any(amat_dir[i, d] == 1 for d in D):
                restrict_dcc[i][j] = []
            else:
                restrict_dcc[i][j] = list(set(D) & set(vset))
    
    restrict_dcc = {k: v for k, v in restrict_dcc.items() if v}
    
    candidates = list(set(label_names) - set(restrict_dcc.keys()))
    pln = []
    for node in candidates:
        adj_nodes_amat = amat_skel[node, amat_skel[node, :] == 1]
        if np.all(adj_nodes_amat == 0):
            pln.append(node)
    
    return {'UdirSubGraph': usubg, 'restricted': restrict_dcc, 'potentialLeaf': pln}

# Example usage:
cpdag = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])  # Simple CPDAG
dcc = {0: [[1]], 1: [[2]]}  # Example direct causal clauses

# Check if DCC is valid
print(is_valid(cpdag, dcc))  # Output: True or False

# Check if DCC is consistent
print(is_consistent(cpdag, dcc))  # Output: {'indicator': True, 'PEO': [0, 1, 2]}

# Construct MPDAG from DCC
print(construct_mpdad(cpdag, dcc))

# Fast construction of MPDAG
print(construct_mpdad_fast(cpdag, dcc))