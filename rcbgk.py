import networkx as nx
import numpy as np
from itertools import combinations
import random
import os
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Set
from graphviz import Digraph
import base64
from pgmpy.base import PDAG


def criticalSet(g, xlb, zlb):
    """Find the critical set of x with respect to z in a chordal graph.
    
    Args:
        g: networkx DiGraph
        xlb: vertex label
        zlb: list of vertex labels
    
    Returns:
        set: Critical set
    """
    nodes = list(g.nodes())
    x = nodes.index(xlb)
    z = [nodes.index(z) for z in zlb]
    n = len(nodes)
    print(f"criticalSet. xlb={xlb}, zlb={zlb}")
    S = [(x, -1, -1)]  # (node, prev, crit)
    C = [] # 使用numpy数组替代矩阵
    dlp = 0
    
    while S:
        dlp += 1
        e = S.pop(0)
        # print(f"criticalSet. x={x}, S={S}, C={C}")
        if e[0] not in z and e[1] != -1:
            # Get neighbors excluding previous node
            neighbors = list(g.neighbors(nodes[e[0]]))
            prev_node = nodes[e[1]]
            if prev_node in neighbors:
                neighbors.remove(prev_node)
            
            # 检查每个邻居
            for alpha in neighbors:
                alpha_idx = nodes.index(alpha)
                # 检查alpha是否与前一个节点相邻
                if prev_node not in g.neighbors(alpha):
                    if e[1] == x:
                        S.append((alpha_idx, e[0], e[0]))
                    else:
                        S.append((alpha_idx, e[0], e[2]))
                        
        elif e[0] not in z:
            # 对于起始节点或其他非z节点
            for alpha in g.neighbors(nodes[e[0]]):
                alpha_idx = nodes.index(alpha)
                S.append((alpha_idx, e[0], alpha_idx))
        else:
            # 当前节点是z中的节点
            C.append(e[2])
            
        if dlp > (n * n + 1):
            break
            
    # 获取非零且唯一的临界节点
    C = list(set(C))
    return {nodes[c] for c in C}

def find_ancestors_and_chaincomp(amat, x):
    """Find ancestors and chain components for a given node in a CPDAG.
    
    Args:
        amat: adjacency matrix of CPDAG (numpy array)
        x: node label name or index
        
    Returns:
        dict with ancestors and chain components
    """
    n = amat.shape[0]
    label_names = list(range(n))
    
    # Convert label to index if needed
    if isinstance(x, str):
        try:
            x = int(x)
        except:
            raise ValueError(f"{x} is not a valid node label")
    
    # Create skeleton and directed/undirected subgraphs
    amat_skel = amat + amat.T
    amat_skel[amat_skel != 0] = 1
    amat_dir = amat - amat.T 
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir

    # Find ancestors using breadth-first search
    w_an = [x]  # Working set for ancestors
    res_an = []  # Result set for ancestors
    while w_an:
        node = w_an.pop(0)
        # Find nodes with directed edges to current node
        an_tmp = np.where(amat_dir[:, node] == 1)[0]
        # Add new ancestors to working set
        w_an.extend([i for i in an_tmp if i not in res_an and i not in w_an])
        res_an.append(node)

    # Find chain components using breadth-first search
    w_cp = [x]  # Working set for chain components
    res_cp = []  # Result set for chain components
    while w_cp:
        node = w_cp.pop(0)
        # Find nodes connected by undirected edges
        cp_tmp = np.where(amat_udir[:, node] == 1)[0]
        # Add new chain component nodes to working set
        w_cp.extend([i for i in cp_tmp if i not in res_cp and i not in w_cp])
        res_cp.append(node)

    return {'an': res_an, 'cp': res_cp}

# ======== transform cbgk =========
def transCbgk(cpdag, bgk):
    """将背景知识转换为直接因果子句
    
    Args:
        cpdag: networkx DiGraph
        bgk: 背景知识列表，每项为 (x, y, t) 其中:
             0: 非因果关系
             1: 因果关系
             
    Returns:
        dict: 直接因果子句
    """
    # 获取邻接矩阵
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    n = len(label_names)
    
    # 创建骨架图和有向/无向子图
    amat_skel = amat + amat.T
    amat_skel[amat_skel != 0] = 1
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir
    
    # 初始化DCC
    dcc = {}
    
    if not bgk:
        return dcc
    print(f"transCbgk. bgk={bgk}")
    for x, y, t in bgk:
        # 转换标签为索引
        x_idx = label_names.index(x) if isinstance(x, str) else x
        y_idx = label_names.index(y) if isinstance(y, str) else y
        x_label = label_names[x_idx]
        y_label = label_names[y_idx]
        
        if amat_skel[x_idx, y_idx] == 1:
            # x与y相邻
            if t == 0:  # 非因果关系
                if y_label not in dcc:
                    dcc[y_label] = []
                dcc[y_label].append({x_label})
            else:  # t == 1, 因果关系
                if x_label not in dcc:
                    dcc[x_label] = []
                dcc[x_label].append({y_label})
        else:
            # x与y不相邻
            # 寻找临界集
            cset = criticalSet(cpdag, x_label, [y_label])
            print(f"transCbgk. cset={cset}")
            if cset:
                if t == 0:  # 非因果关系
                    for element in cset:
                        if element not in dcc:
                            dcc[element] = []
                        dcc[element].append({x_label})
                else:  # t == 1, 因果关系
                    if x_label not in dcc:
                        dcc[x_label] = []
                    dcc[x_label].append(cset)
            elif t == 1:  # 直接因果关系但无路径 - 矛盾
                if x_label not in dcc:
                    dcc[x_label] = []
                dcc[x_label].append(set())
    
    return {k:v for k,v in dcc.items() if v}


# ======== random generation ========
def gen_valid_dcc(cpdag, num):
    """
    Generate valid direct causal clauses
    Args:
        cpdag: networkx Graph
        num: number of clauses to generate
    Returns:
        dict of valid direct causal clauses
    """
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    n = len(label_names)
    
    # Create directed/undirected subgraphs
    amat_skel = amat + amat.T
    amat_skel[amat_skel != 0] = 1
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir
    
    # Randomly sample tails
    tails = np.random.choice(label_names, num)
    uniq_tails, counts = np.unique(tails, return_counts=True)
    dcc = {tail: [] for tail in uniq_tails}
    
    for i, x in enumerate(uniq_tails):
        x_idx = label_names.index(x)
        adj_x = np.where(amat_skel[x_idx] == 1)[0]
        pa_x = np.where(amat_dir[:, x_idx] == 1)[0]
        ch_x = np.where(amat_dir[x_idx] == 1)[0]
        
        if all(a in pa_x for a in adj_x):
            continue
            
        x_num = counts[i]
        for _ in range(x_num):
            # Randomly select adjacent nodes
            adj_x_labels = [label_names[j] for j in adj_x]
            ind = np.random.randint(2, size=len(adj_x_labels))
            select = [adj_x_labels[j] for j in range(len(adj_x_labels)) if ind[j]]
            
            new_select = []
            if all(s in [label_names[j] for j in pa_x] for s in select):
                can = [label_names[j] for j in adj_x if j not in pa_x]
                while not new_select:
                    new_select = [c for c in can if np.random.random() < 0.5]
                    
            if not all(s in [label_names[j] for j in ch_x] for s in select):
                if np.random.random() < 0.8:
                    select = [s for s in select if s not in [label_names[j] for j in ch_x]]
                    
            D = list(set(select + new_select))
            if D:  # Only append non-empty lists
                dcc[x].append(D)
                
        # Remove duplicates
        if dcc[x]:
            dcc[x] = [list(x) for x in set(tuple(sorted(x)) for x in dcc[x])]
            
    return {k:v for k,v in dcc.items() if v}

def gen_invalid_dcc(cpdag, num, prob=0.9):
    """
    Generate invalid direct causal clauses
    Args:
        cpdag: networkx Graph
        num: number of clauses to generate
        prob: probability of adding invalid nodes
    Returns:
        dict of invalid direct causal clauses
    """
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    n = len(label_names)
    
    # Create directed/undirected subgraphs
    amat_skel = amat + amat.T
    amat_skel[amat_skel != 0] = 1
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir
    
    # Randomly sample tails
    tails = np.random.choice(label_names, num)
    uniq_tails, counts = np.unique(tails, return_counts=True)
    dcc = {tail: [] for tail in uniq_tails}
    
    for i, x in enumerate(uniq_tails):
        x_idx = label_names.index(x)
        adj_x = np.where(amat_skel[x_idx] == 1)[0]
        pa_x = np.where(amat_dir[:, x_idx] == 1)[0]
        ch_x = np.where(amat_dir[x_idx] == 1)[0]
        
        x_num = counts[i]
        for _ in range(x_num):
            # Randomly select adjacent nodes
            adj_x_labels = [label_names[j] for j in adj_x]
            ind = np.random.randint(2, size=len(adj_x_labels))
            select = [adj_x_labels[j] for j in range(len(adj_x_labels)) if ind[j]]
            
            new_select = []
            if np.random.random() < prob:
                new_select = [np.random.choice(label_names)]
                
            D = list(set(select + new_select))
            if D:  # Only append non-empty lists
                dcc[x].append(D)
                
        # Remove duplicates
        if dcc[x]:
            dcc[x] = [list(x) for x in set(tuple(sorted(x)) for x in dcc[x])]
            
    return {k:v for k,v in dcc.items() if v}

# ======== validity and consistency ========
def is_valid(cpdag, dcc):
    """
    Check if direct causal clauses are valid
    Args:
        cpdag: networkx Graph
        dcc: dict of direct causal clauses
    Returns:
        bool indicating validity
    """
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    
    # Create directed/undirected subgraphs
    amat_skel = amat + amat.T
    amat_skel[amat_skel != 0] = 1
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir
    
    if dcc:
        for x, clauses in dcc.items():
            x_idx = label_names.index(x)
            for D in clauses:
                if not D:  # Empty clause
                    return False
                    
                D_idx = [label_names.index(d) for d in D]
                # Check if all nodes in D are adjacent to x
                if any(amat_skel[x_idx, d] == 0 for d in D_idx):
                    return False
                    
                # Check if all nodes in D are parents of x
                if all(amat_dir[d, x_idx] == 1 for d in D_idx):
                    return False
                    
    return True

def is_consistent(cpdag, dcc):
    """
    Check if direct causal clauses are consistent
    Args:
        cpdag: networkx Graph
        dcc: dict of direct causal clauses
    Returns:
        dict with consistency indicator and perfect elimination ordering
    """
    if not is_valid(cpdag, dcc):
        return {'indicator': False, 'PEO': None}
        
    label_names = list(cpdag.nodes())
    vset = label_names.copy()
    peo = []
    
    while vset:
        pln = restriction(cpdag, vset, dcc)['potential_leaf']
        if not pln:
            return {'indicator': False, 'PEO': None}
        else:
            vset.remove(pln[0])
            peo.append(pln[0])
            
    return {'indicator': True, 'PEO': peo}

def restriction(cpdag, vset, dcc):
    """Get restriction of graph and DCC to vertex subset.
    
    Args:
        cpdag: networkx Graph
        vset: list of vertex labels
        dcc: dict of direct causal clauses
        
    Returns:
        dict with restricted graph and DCC
    """
    if not vset:
        raise ValueError("vset should be non-empty.")
        
    # Get adjacency matrices
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    
    # Create directed/undirected subgraphs
    amat_skel = amat + amat.T
    amat_skel[amat_skel != 0] = 1
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir
    
    # Get indices for vset
    vset_idx = [label_names.index(v) for v in vset]
    
    # Create restricted undirected subgraph
    amat_udir_sub = amat_udir[np.ix_(vset_idx, vset_idx)]
    g_sub = nx.from_numpy_array(amat_udir_sub)
    nx.relabel_nodes(g_sub, {i: vset[i] for i in range(len(vset))}, copy=False)
    
    # Restrict DCC
    restricted_dcc = {k: v.copy() for k, v in dcc.items() if k in vset}
    
    if restricted_dcc:
        for x, clauses in list(restricted_dcc.items()):
            x_idx = label_names.index(x)
            new_clauses = []
            
            for D in clauses:
                # Check if x has directed edges to any node in D
                if any(amat_dir[x_idx, label_names.index(d)] == 1 for d in D):
                    continue
                    
                # Restrict x ⊥ D to G*_u
                D_restricted = {d for d in D 
                              if d in vset and 
                              amat_udir[x_idx, label_names.index(d)] == 1}
                
                if D_restricted:
                    new_clauses.append(D_restricted)
                    
            if new_clauses:
                restricted_dcc[x] = new_clauses
            else:
                del restricted_dcc[x]
                
    # Find potential leaf nodes
    candidates = [v for v in vset if v not in restricted_dcc]
    pln = []
    
    if candidates:
        for c in candidates:
            c_idx = vset.index(c)
            adj_nodes = np.where(amat_udir_sub[c_idx] == 1)[0]
            
            # Check if node forms a complete subgraph with its neighbors
            if len(adj_nodes) <= 1:
                pln.append(c)
            else:
                adj_submat = amat_udir_sub[np.ix_(adj_nodes, adj_nodes)]
                # 修改判断条件：检查是否所有元素都为0（除了对角线）
                if np.sum(adj_submat == 0) == len(adj_nodes):
                    pln.append(c)
                    
    return {
        'udir_subgraph': g_sub,
        'restricted': restricted_dcc,
        'potential_leaf': pln
    }

# ======== construct MPDAG ========
def construct_moc(cpdag, x, dcc, simplified=True):
    """构造节点x的最大可定向组件
    
    Args:
        cpdag: networkx Graph
        x: 顶点标签
        dcc: 直接因果子句字典
        simplified: 如果为True，仅返回相邻节点
    Returns:
        MOC中的节点列表或包含MOC和相邻节点的字典
    """
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    n = len(label_names)
    
    # 创建有向/无向子图
    amat_skel = amat + amat.T
    amat_skel[amat_skel != 0] = 1
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir
    print(f"construct_moc. x={x}, dcc={dcc}")
    if not is_valid(cpdag, dcc):
        raise ValueError("The direct causal clauses are not valid with respect to the CPDAG.")
        
    vset = label_names.copy()
    while vset:
        rst = restriction(cpdag, vset, dcc)
        pln = rst['potential_leaf']
        print(f"vset={vset}, pln={pln}")
        if not pln:
            raise ValueError("The direct causal clauses are inconsistent with respect to the CPDAG.")
        elif len(pln) == 1 and pln[0] == x:
            moc_set = vset
            break
        else:
            # 关键修改：使用setdiff而不是next
            to_remove = list(set(pln) - {x})[0]
            vset.remove(to_remove)
    print(f"vset={vset}")
    # 获取无向邻居
    adj = [v for v in vset if amat_udir[label_names.index(x), label_names.index(v)] == 1]
    
    return adj if simplified else {'MOC': moc_set, 'adj': adj}

def construct_mpdag(cpdag, dcc):
    """从CPDAG和直接因果子句构造MPDAG
    
    Args:
        cpdag: networkx DiGraph
        dcc: 直接因果子句字典
        
    Returns:
        networkx DiGraph: 构造的MPDAG
    """
    if not dcc:
        return cpdag.copy()
        
    # 获取邻接矩阵
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    n = len(label_names)
    print(f"amat: {amat}")
    print(f"label_names: {label_names}")
    # 创建有向/无向子图
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir
    
    # 处理每个DCC
    for x in label_names:
        moc = construct_moc(cpdag, x, dcc)
        print(f"x: {x}, moc: {moc}")
        if moc:  # 如果找到MOC
            x_idx = label_names.index(x)
            for y in moc:
                y_idx = label_names.index(y)
                # 将边定向为 y -> x
                amat[y_idx, x_idx] = 1
                amat[x_idx, y_idx] = 0
    
    
    mpdag = nx.DiGraph()
    mpdag.add_nodes_from(label_names)
    for i in range(n):
        for j in range(n):
            if amat[i,j] == 1:
                mpdag.add_edge(label_names[i], label_names[j])
    
    return mpdag

def is_valid_dcc(cpdag, dcc):
    """Check if direct causal clauses are valid.
    
    Args:
        cpdag: A CPDAG as networkx DiGraph
        dcc: Dictionary of direct causal clauses
        
    Returns:
        True if valid, False otherwise
    """
    if not dcc:
        return True
        
    # Get adjacency matrix
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    
    # Create directed/undirected subgraphs
    amat_skel = amat + amat.T
    amat_skel[amat_skel != 0] = 1
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    
    for x, clauses in dcc.items():
        x_idx = label_names.index(x)
        for D in clauses:
            # Empty clause is invalid
            if not D:
                return False
                
            # Check adjacency
            for d in D:
                d_idx = label_names.index(d)
                if amat_skel[x_idx, d_idx] == 0:
                    return False
            
            # Check if all vertices in D are already parents
            D_idx = [label_names.index(d) for d in D]
            if all(amat_dir[d, x_idx] == 1 for d in D_idx):
                return False
    
    return True

# ======== IDA ========
def bgk_ida(x, y, mcov, cpdag, mpdag=None, bgk=None, verbose=False):
    """IDA algorithm with background knowledge.
    
    Args:
        x, y: vertex labels or indices
        mcov: covariance matrix
        cpdag: networkx DiGraph
        mpdag: optional pre-computed MPDAG
        bgk: background knowledge matrix
        verbose: print details
        
    Returns:
        array of possible causal effects
    """
    label_names = list(cpdag.nodes())
    
    # Convert labels to indices if needed
    if isinstance(x, str):
        x = label_names.index(x)
    if isinstance(y, str):
        y = label_names.index(y)
        
    # Check if CPDAG is valid
    amat_cpdag = nx.to_numpy_array(cpdag)
    amat_cpdag[amat_cpdag != 0] = 1
    
    if bgk is None:
        beta_hat = ida(x, y, mcov, cpdag)
    else:
        # Transform background knowledge to DCC
        dcc = transCbgk(cpdag, bgk)
        
        # Construct MPDAG if not provided
        if mpdag is None:
            mpdag = construct_mpdag(cpdag, dcc)
            
        # Get adjacency matrix
        amat = nx.to_numpy_array(mpdag)
        
        if verbose:
            print(f"MPDAG:\n{amat}")
            print(f"DCC: {dcc}")
            
        # Get definite parents and siblings
        wgt_est = (amat != 0)
        wgt_unique = np.logical_and(wgt_est, ~wgt_est.T)
        pa1 = np.where(wgt_unique[:,x])[0]
        
        if y in pa1:
            beta_hat = [0]
        else:
            wgt_ambig = np.logical_and(wgt_est, wgt_est.T)
            pa2 = np.where(wgt_ambig[x,:])[0]
            
            if verbose:
                print(f"\nx={x} y={y}\npa1={pa1}\npa2={pa2}")
                
            if len(pa2) == 0:
                beta_hat = [lm_cov(mcov, y, np.append(x, pa1))]
                if verbose:
                    print(f"Fit - y:{y} x:{np.append(x, pa1)} |b_hat={beta_hat[0]}")
                    
            else:
                beta_hat = []
                
                # Try empty pa2_t
                pa2_f = pa2
                pa2_t = np.array([], dtype=int)
                if is_locally_consistent(cpdag, dcc, label_names[x], 
                                          [label_names[i] for i in pa2_t],
                                          [label_names[i] for i in pa2_f]):
                    b = lm_cov(mcov, y, np.append(x, pa1))
                    beta_hat.append(b)
                    if verbose:
                        print(f"Fit - y:{y} x:{np.append(x, pa1)} |b_hat={b}")
                        
                # Try single nodes in pa2_t
                for i2 in range(len(pa2)):
                    pa2_f = np.delete(pa2, i2)
                    pa2_t = np.array([pa2[i2]])
                    
                    if is_locally_consistent(cpdag, dcc, label_names[x],
                                          [label_names[i] for i in pa2_t],
                                          [label_names[i] for i in pa2_f]):
                        if y in pa2_t:
                            beta_hat.append(0)
                        else:
                            b = lm_cov(mcov, y, np.append(np.append(x, pa1), pa2[i2]))
                            beta_hat.append(b)
                            if verbose:
                                print(f"Fit - y:{y} x:{np.append(np.append(x, pa1), pa2[i2])} |b_hat={b}")
                                
                # Try larger subsets of pa2_t
                if len(pa2) > 1:
                    for i in range(2, len(pa2) + 1):
                        for pa_tmp in combinations(pa2, i):
                            pa2_t = np.array(pa_tmp)
                            pa2_f = np.array([p for p in pa2 if p not in pa2_t])
                            
                            if is_locally_consistent(cpdag, dcc, label_names[x],
                                                  [label_names[i] for i in pa2_t],
                                                  [label_names[i] for i in pa2_f]):
                                if y in pa2_t:
                                    beta_hat.append(0)
                                else:
                                    b = lm_cov(mcov, y, np.append(np.append(x, pa1), pa2_t))
                                    beta_hat.append(b)
                                    if verbose:
                                        print(f"Fit - y:{y} x:{np.append(np.append(x, pa1), pa2_t)} |b_hat={b}")
                                        
    return np.array(beta_hat)

def has_extension(amat, amat_skel, x, pa1, pa2_t, pa2_f):
    """Check if the orientation creates a valid extension.
    
    Args:
        amat: adjacency matrix
        amat_skel: skeleton matrix
        x: node index
        pa1: definite parents
        pa2_t: potential parents
        pa2_f: potential children
    
    Returns:
        bool: True if orientation is valid
    """
    # Check if undirected edges that are pointed to x create a new v-structure
    # or directed triangle containing x
    res = False
    
    # Check v-structure
    if len(pa2_t) > 0:
        # Check whether all pa1 and pa2_t are connected
        if len(pa1) > 0:
            res = min([amat_skel[p1, p2] for p1 in pa1 for p2 in pa2_t]) == 0
        
        # All pa2_t have to be connected
        if not res and len(pa2_t) > 1:
            A2 = np.array([[amat_skel[i,j] for j in pa2_t] for i in pa2_t])
            np.fill_diagonal(A2, 1)
            res = min(A2.flatten()) == 0
            
    # 由于 Meek R1，这部分检查可以移除
    # if not res and len(pa2_f) > 0:
    #     A = amat - amat.T
    #     A[A < 0] = 0
    #     papa = set()
    #     for node in pa2_f:
    #         papa.update(np.where(A[:,node] != 0)[0])
    #     papa = papa - {x}
    #     if len(papa) > 0:
    #         res = min([amat_skel[x,p] for p in papa]) == 0
            
    # Check directed triangle containing x
    if not res:
        A = amat - amat.T
        A[A < 0] = 0
        nA = A.copy()
        
        # 修正：pa2_t 是潜在父节点，应该指向 x
        for p in pa2_t:
            nA[x,p] = 0
            nA[p,x] = 1  # p -> x
            
        # pa2_f 是潜在子节点，应该被 x 指向
        for c in pa2_f:
            nA[x,c] = 1  # x -> c
            nA[c,x] = 0
            
        cd = np.where(nA[:,x] != 0)[0]  # 指向 x 的节点
        pa = np.where(nA[x,:] != 0)[0]  # x 指向的节点
        
        if len(cd) > 0 and len(pa) > 0:
            subA = nA[np.ix_(pa,cd)]
            res = max(subA.flatten()) == 1
            
    return not res

def ida(x, y, mcov, cpdag, method="local", verbose=False):
    """Standard IDA algorithm without background knowledge.
    
    Args:
        x: vertex index
        y: vertex index
        mcov: covariance matrix
        cpdag: networkx Graph
        method: "local" or "global"
        verbose: print details if True
        
    Returns:
        array of possible causal effects
    """
    # Get adjacency matrix
    amat = nx.to_numpy_array(cpdag)
    wgt_est = (amat != 0).astype(int)
    
    # Get directed and undirected edges
    tmp = wgt_est - wgt_est.T
    tmp[tmp < 0] = 0
    wgt_unique = tmp
    
    # Get definite parents (已修正)
    pa1 = np.where(wgt_unique[:,x])[0]
    
    if y in pa1:
        return np.array([0])
        
    # Get potential parents/children
    wgt_ambig = wgt_est - wgt_unique
    pa2 = np.where(wgt_ambig[x,:])[0]  # 修改这里：应该看行，因为是无向边
    
    if verbose:
        print(f"\nx={x} y={y}\npa1={pa1}\npa2={pa2}\n")
        
    if len(pa2) == 0:
        beta_hat = lm_cov(mcov, y, np.append(x, pa1))
        if verbose:
            print(f"Fit - y:{y} x:{np.append(x, pa1)} |b_hat={beta_hat}")
        return np.array([beta_hat])
        
    # Get skeleton
    amat_skel = amat + amat.T
    amat_skel[amat_skel != 0] = 1
    
    beta_hat = []
    
    # Try empty pa2_t
    pa2_t = np.array([])
    pa2_f = pa2
    if has_extension(amat, amat_skel, x, pa1, pa2_t, pa2_f):
        b = lm_cov(mcov, y, np.append(x, pa1))
        beta_hat.append(b)
        if verbose:
            print(f"Fit - y:{y} x:{np.append(x, pa1)} |b_hat={b}")
            
    # Try single nodes in pa2_t
    for i2 in range(len(pa2)):
        pa2_f = np.delete(pa2, i2)
        pa2_t = np.array([pa2[i2]])
        
        if has_extension(amat, amat_skel, x, pa1, pa2_t, pa2_f):
            if y in pa2_t:
                beta_hat.append(0)
            else:
                b = lm_cov(mcov, y, np.append(np.append(x, pa1), pa2[i2]))
                beta_hat.append(b)
                if verbose:
                    print(f"Fit - y:{y} x:{np.append(np.append(x, pa1), pa2[i2])} |b_hat={b}")
                    
    # Try larger subsets of pa2_t
    if len(pa2) > 1:
        for i in range(2, len(pa2) + 1):
            for pa_tmp in combinations(pa2, i):
                pa2_t = np.array(pa_tmp)
                pa2_f = np.array([p for p in pa2 if p not in pa2_t])
                
                if has_extension(amat, amat_skel, x, pa1, pa2_t, pa2_f):
                    if y in pa2_t:
                        beta_hat.append(0)
                    else:
                        b = lm_cov(mcov, y, np.append(np.append(x, pa1), pa2_t))
                        beta_hat.append(b)
                        if verbose:
                            print(f"Fit - y:{y} x:{np.append(np.append(x, pa1), pa2_t)} |b_hat={b}")
                            
    return np.array(beta_hat)

def lm_cov(C, y, x):
    """
    Calculate regression coefficient using covariance matrix
    Args:
        C: covariance matrix
        y: response variable index
        x: predictor variable indices
    Returns:
        float regression coefficient
    """
    return float(np.linalg.solve(C[np.ix_(x, x)], C[np.ix_(x, [y])])[0])

def is_globally_consistent(cpdag, dcc, x, pa2_t, pa2_f):
    """Check if orientation is globally consistent.
    
    Args:
        cpdag: networkx Graph representing CPDAG
        dcc: dict of direct causal clauses
        x: vertex label
        pa2_t: list of vertex labels (parents)
        pa2_f: list of vertex labels (children)
        
    Returns:
        bool indicating consistency
    """
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    
    # Create directed/undirected subgraphs
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir
    
    # Get vertices in undirected component
    vset = [x]
    for v in label_names:
        if amat_udir[label_names.index(x), label_names.index(v)] == 1:
            vset.append(v)
    
    # Create new DCC with orientations
    dcc_new = dcc.copy()
    
    # Add orientations towards x
    for p in pa2_t:
        if p not in dcc_new:
            dcc_new[p] = []
        dcc_new[p].append({x})
        
    # Add orientations away from x
    for c in pa2_f:
        if x not in dcc_new:
            dcc_new[x] = []
        dcc_new[x].append({c})
    
    # Get restriction of graph and DCC
    rst = restriction(cpdag, vset, dcc_new)
    
    # Check consistency
    return is_consistent(rst['udir_subgraph'], rst['restricted'])['indicator']

def effect_distance(x, y, tdag, effect_pair):
    """
    Compute distance between true and estimated causal effects
    Args:
        x, y: vertex labels or indices
        tdag: true DAG (networkx DiGraph)
        effect_pair: array of estimated effects
    Returns:
        tuple of (min_dist, max_dist, mse)
    """
    # Get true adjacency matrix and covariance
    true_mat = nx.to_numpy_array(tdag)
    mcov = true_cov(tdag)
    label_names = list(tdag.nodes())
    
    # Convert labels to indices if needed
    if isinstance(x, str):
        x_idx = label_names.index(x)
    else:
        x_idx = x
        
    if isinstance(y, str):
        y_idx = label_names.index(y)
    else:
        y_idx = y
    
    # Get parents and descendants
    pa = list(tdag.predecessors(label_names[x_idx]))
    pa_idx = [label_names.index(p) for p in pa]
    de = nx.descendants(tdag, label_names[x_idx])
    
    # Calculate true effect
    if label_names[y_idx] in de:
        true_eff = lm_cov(mcov, y_idx, np.append(x_idx, pa_idx))
    else:
        true_eff = 0
        
    # Calculate distances
    dists = np.abs(effect_pair - true_eff)
    mse = np.mean(dists**2)
    
    return (np.min(dists), np.max(dists), mse)

def int_mse(x, y, tdag, effect_pair):
    """
    Calculate intervention MSE
    Args:
        x, y: vertex labels or indices
        tdag: true DAG (networkx DiGraph)
        effect_pair: array of estimated effects
    Returns:
        intervention MSE
    """
    # Get true adjacency matrix and covariance
    true_mat = nx.to_numpy_array(tdag)
    mcov = true_cov(tdag)
    label_names = list(tdag.nodes())
    
    # Convert indices to labels if needed
    if isinstance(x, int):
        x_idx = x
    else:
        x_idx = label_names.index(x)
        
    if isinstance(y, int):
        y_idx = y
    else:
        y_idx = label_names.index(y)
    
    # Get parents and descendants
    pa = np.where(true_mat[:, x_idx] != 0)[0]
    de = nx.descendants(tdag, x_idx)
    
    # Calculate true effect
    if y_idx in de:
        true_eff = lm_cov(mcov, y_idx, np.append(x_idx, pa))
    else:
        true_eff = 0
        
    # Calculate intervention MSE
    ub = np.max(effect_pair)
    lb = np.min(effect_pair)
    
    if true_eff > ub:
        return true_eff - ub
    elif true_eff < lb:
        return lb - true_eff
    else:
        return 0

def bgk_ida_semilocal(x, y, mcov, dag, bgk=None, verbose=False):
    """
    Calculate the possible causal effects with background knowledge
    Args:
        x: source node
        y: target node
        mcov: covariance matrix
        dag: networkx DiGraph
        bgk: background knowledge (optional)
        verbose: print details (default: False)
    Returns:
        numpy array of possible causal effects
    """
    # Get CPDAG and DCC
    cpdag = dag_to_cpdag(dag)
    dcc = bgk_to_dcc(cpdag, bgk) if bgk is not None else {}
    
    # Get label names and convert to indices if needed
    label_names = list(dag.nodes())
    if isinstance(x, str):
        x = label_names.index(x)
    if isinstance(y, str):
        y = label_names.index(y)
    x_lb = label_names[x]
    
    # Get adjacency matrix
    amat = nx.to_numpy_array(cpdag)
    
    if verbose:
        print(f"CPDAG:\n{amat}")
        print(f"DCC: {dcc}")
        
    # Get definite parents and siblings
    wgt_est = (amat != 0)
    wgt_unique = np.logical_and(wgt_est, ~wgt_est.T)
    pa1 = np.where(wgt_unique[:,x])[0]
    
    if y in pa1:
        beta_hat = [0]
    else:
        wgt_ambig = np.logical_and(wgt_est, wgt_est.T)
        pa2 = np.where(wgt_ambig[x,:])[0]
        
        if verbose:
            print(f"\nx={x} y={y}\npa1={pa1}\npa2={pa2}")
            
        if len(pa2) == 0:
            beta_hat = [lm_cov(mcov, y, np.append(x, pa1))]
            if verbose:
                print(f"Fit - y:{y} x:{np.append(x, pa1)} |b_hat={beta_hat[0]}")
                    
        else:
            beta_hat = []
            
            # Try empty pa2_t
            pa2_f = pa2
            pa2_t = np.array([], dtype=int)
            if is_globally_consistent(cpdag, dcc, x_lb, pa2_t, [label_names[i] for i in pa2_f]):
                b = lm_cov(mcov, y, np.append(x, pa1))
                beta_hat.append(b)
                if verbose:
                    print(f"Fit - y:{y} x:{np.append(x, pa1)} |b_hat={b}")
                    
            # Try single nodes in pa2_t
            for i2 in range(len(pa2)):
                pa2_f = np.delete(pa2, i2)
                pa2_t = np.array([pa2[i2]])
                
                if is_globally_consistent(cpdag, dcc, x_lb, 
                                           [label_names[i] for i in pa2_t],
                                           [label_names[i] for i in pa2_f]):
                    if y in pa2_t:
                        beta_hat.append(0)
                    else:
                        b = lm_cov(mcov, y, np.append(np.append(x, pa1), pa2[i2]))
                        beta_hat.append(b)
                        if verbose:
                            print(f"Fit - y:{y} x:{np.append(np.append(x, pa1), pa2[i2])} |b_hat={b}")
                                
            # Try larger subsets of pa2_t
            if len(pa2) > 1:
                for i in range(2, len(pa2) + 1):
                    for pa_tmp in combinations(pa2, i):
                        pa2_t = np.array(pa_tmp)
                        pa2_f = np.array([p for p in pa2 if p not in pa2_t])
                        
                        if is_globally_consistent(cpdag, dcc, x_lb,
                                                [label_names[i] for i in pa2_t],
                                                [label_names[i] for i in pa2_f]):
                            if y in pa2_t:
                                beta_hat.append(0)
                            else:
                                b = lm_cov(mcov, y, np.append(np.append(x, pa1), pa2_t))
                                beta_hat.append(b)
                                if verbose:
                                    print(f"Fit - y:{y} x:{np.append(np.append(x, pa1), pa2_t)} |b_hat={b}")
                                        
    return np.array(beta_hat)

def true_cov(dag, weights=None):
    """Calculate the covariance matrix of a linear Gaussian DAG model.
    
    Args:
        dag: networkx DiGraph or adjacency matrix
        weights: Optional edge weights matrix
        
    Returns:
        numpy array representing the covariance matrix
    """
    # Convert input to adjacency matrix if needed
    if isinstance(dag, nx.DiGraph):
        B = nx.to_numpy_array(dag)
    else:
        B = np.array(dag)
    
    n = B.shape[0]
    
    # Apply weights if provided
    if weights is not None:
        B = B * weights
    
    # Calculate covariance matrix using the formula Σ = (I-B)^(-1) (I-B)^(-T)
    try:
        IminusB = np.eye(n) - B.T
        IminusB_inv = np.linalg.inv(IminusB)
        cov = IminusB_inv @ IminusB_inv.T
        return cov
    except np.linalg.LinAlgError:
        print("Warning: Matrix is not invertible, returning identity matrix")
        return np.eye(n)

def is_locally_consistent(cpdag, dcc, x, pa2_t, pa2_f):
    """
    Check if orientation is locally consistent
    Args:
        cpdag: networkx Graph
        dcc: dict of direct causal clauses
        x: vertex label
        pa2_t: list of vertex labels (parents)
        pa2_f: list of vertex labels (children)
    Returns:
        bool indicating consistency
    """
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    
    # Create directed/undirected subgraphs
    amat_dir = amat - amat.T
    amat_dir[amat_dir < 0] = 0
    amat_udir = amat - amat_dir
    
    # Get vertices in undirected component
    vset = [x]
    for v in label_names:
        if amat_udir[label_names.index(x), label_names.index(v)] == 1:
            vset.append(v)
    
    # Check local consistency first
    for p in pa2_t:
        if not any(amat_udir[label_names.index(p), label_names.index(v)] == 1 for v in vset):
            return False
            
    for c in pa2_f:
        if not any(amat_udir[label_names.index(c), label_names.index(v)] == 1 for v in vset):
            return False
    
    # Create new DCC with orientations
    dcc_new = dcc.copy()
    
    # Add orientations
    for p in pa2_t:
        if p not in dcc_new:
            dcc_new[p] = []
        dcc_new[p].append({x})
        
    for c in pa2_f:
        if x not in dcc_new:
            dcc_new[x] = []
        dcc_new[x].append({c})
    
    # Check global consistency
    rst = restriction(cpdag, vset, dcc_new)
    return is_consistent(rst['udir_subgraph'], rst['restricted'])['indicator']

def tporder(edges, p):
    """
    Compute topological ordering of nodes in DAG
    Args:
        edges: list of [from, to] edges
        p: number of nodes
    Returns:
        list of nodes in topological order
    """
    from collections import deque
    
    # Initialize data structures
    des_dag = [[] for _ in range(p)]
    par_num = [0 for _ in range(p)]
    
    # Count parents and build descendant lists
    for edge in edges:
        des_dag[edge[0]].append(edge[1])
        par_num[edge[1]] += 1
        
    # Find nodes with zero in-degree
    zero_indeg = deque([i for i in range(p) if par_num[i] == 0])
    
    # Build topological ordering
    s = 0
    tp_ord = [0 for _ in range(p)]
    while zero_indeg:
        v = zero_indeg.popleft()
        tp_ord[s] = v
        s += 1
        children = des_dag[v]
        if not children:
            continue
        for j in children:
            par_num[j] -= 1
            if par_num[j] == 0:
                zero_indeg.appendleft(j)
                
    return tp_ord

def edge_order(edges, p):
    """
    Order edges according to Chickering's algorithm
    Args:
        edges: list of [from, to] edges
        p: number of nodes
    Returns:
        ordered list of edges
    """
    # Get topological ordering
    tp_ord = tporder(edges, p)
    
    # Create mapping from node to position in topological order
    nord_map = [0 for _ in range(p)]
    for i in range(p):
        nord_map[tp_ord[i]] = i
        
    # Sort edges by (child position, -parent position)
    return sorted(edges, key=lambda x: (nord_map[x[1]], -nord_map[x[0]]))

def label_edges(edges, p):
    """
    Label edges as compelled (1) or reversible (-1)
    Args:
        edges: list of [from, to] edges
        p: number of nodes
    Returns:
        list of [from, to, label] edges
    """
    # Initialize data structures
    ans_dag = [[] for _ in range(p)]
    edge_num = [[] for _ in range(p)]
    edge_ord = edge_order(edges, p)
    
    # Build parent lists and edge indices
    for i, edge in enumerate(edge_ord):
        ans_dag[edge[1]].append(edge[0])
        edge_num[edge[1]].append(i)
        
    # Label edges
    unlabeled = len(edges)
    edge_label = [0 for _ in range(unlabeled)]
    
    while unlabeled > 0:
        v = edge_label.index(0)
        x = edge_ord[v][0]
        y = edge_ord[v][1]
        
        # Check both parents and children
        is_compelled = False
        
        # Check parents of x that are compelled
        for wi in [i for i in range(len(ans_dag[x])) if edge_label[edge_num[x][i]] == 1]:
            w = ans_dag[x][wi]
            if w not in ans_dag[y]:
                is_compelled = True
                break
            else:
                # Label edge from w to y as compelled
                wyi = ans_dag[y].index(w)
                if edge_label[edge_num[y][wyi]] == 0:
                    unlabeled -= 1
                edge_label[edge_num[y][wyi]] = 1
        
        if is_compelled:
            # Label all unlabeled edges into y as compelled
            for j in edge_num[y]:
                if edge_label[j] == 0:
                    unlabeled -= 1
                edge_label[j] = 1
            continue
        
        # Check if y has any parents that are not parents of x
        if any(z != x and z not in ans_dag[x] for z in ans_dag[y]):
            # Label all unlabeled edges into y as compelled
            for yi in edge_num[y]:
                if edge_label[yi] == 0:
                    unlabeled -= 1
                    edge_label[yi] = 1
        else:
            # Label all unlabeled edges into y as reversible
            for yi in edge_num[y]:
                if edge_label[yi] == 0:
                    unlabeled -= 1
                    edge_label[yi] = -1
                    
    # Add labels to edges
    for i in range(len(edges)):
        edge_ord[i].append(edge_label[i])
        
    return edge_ord

def dag_to_cpdag(dag):
    """
    Convert a DAG to its corresponding CPDAG using Chickering's algorithm
    Args:
        dag: networkx DiGraph
    Returns:
        networkx Graph representing the CPDAG
    """
    # Get edges from DAG and node labels
    label_names = list(dag.nodes())
    edges = list(dag.edges())
    edges = [[label_names.index(u), label_names.index(v)] for u, v in edges]
    p = len(label_names)
    
    # Label edges as compelled or reversible
    labeled_edges = label_edges(edges, p)
    
    # Create CPDAG
    cpdag = nx.DiGraph()
    cpdag.add_nodes_from(label_names)
    
    # Add edges to CPDAG
    for edge in labeled_edges:
        u, v, label = edge
        u_name = label_names[u]
        v_name = label_names[v]
        if label == 1:  # Compelled edge
            cpdag.add_edge(u_name, v_name)
        else:  # Reversible edge
            cpdag.add_edge(u_name, v_name)
            cpdag.add_edge(v_name, u_name)
            
    return cpdag

def generate_undirected_chordal(n, edge_ratio):
    """Generate an undirected chordal graph with n nodes and specified edge ratio.
    Returns a networkx Graph object to be consistent with other functions in rcbgk.
    
    Parameters:
    -----------
    n : int
        Number of nodes
    edge_ratio : float
        Ratio of edges to maximum possible edges (between 0 and 1)
        
    Returns:
    --------
    nx.Graph
        An undirected chordal graph
    """
    # Calculate number of edges
    max_edges = n * (n - 1) // 2
    num_edges = min(int(max_edges * edge_ratio), max_edges)
    
    # Start with a random tree
    G = nx.random_tree(n)
    
    # Add edges while maintaining chordal property
    edges_to_add = num_edges - (n - 1)  # Tree already has n-1 edges
    potential_edges = [(i, j) for i in range(n) for j in range(i+1, n) 
                      if not G.has_edge(i, j)]
    
    while edges_to_add > 0 and potential_edges:
        # Randomly select an edge
        i, j = random.choice(potential_edges)
        potential_edges.remove((i, j))
        
        # Try adding the edge
        G.add_edge(i, j)
        
        # Check if graph remains chordal
        if nx.is_chordal(G):
            edges_to_add -= 1
        else:
            G.remove_edge(i, j)
    
    # Relabel nodes to match rcbgk convention (using letters)
    mapping = {i: chr(65 + i) for i in range(n)}
    G = nx.relabel_nodes(G, mapping)
    
    return G

def generate_undirected_chordal_set(n, edge_ratios, num_graphs=10):
    """Generate a set of undirected chordal graphs with different edge ratios.
    
    Parameters:
    -----------
    n : int
        Number of nodes
    edge_ratios : list of float
        List of edge ratios to test
    num_graphs : int
        Number of graphs to generate for each ratio
        
    Returns:
    --------
    dict
        Dictionary containing the generated graphs and their properties
    """
    results = []
    
    for ratio in edge_ratios:
        for i in range(num_graphs):
            G = generate_undirected_chordal(n, ratio)
            
            # Calculate graph properties
            result = {
                'graph': G,
                'graph_id': i,
                'edge_ratio': ratio,
                'num_edges': G.number_of_edges(),
                'density': nx.density(G),
                'max_clique_size': len(max(nx.find_cliques(G), key=len)),
                'num_cliques': len(list(nx.find_cliques(G))),
                'is_chordal': nx.is_chordal(G),
                'is_connected': nx.is_connected(G)
            }
            results.append(result)
    
    return results

def analyze_chordal_graphs(graphs_data):
    """Analyze properties of generated chordal graphs.
    
    Parameters:
    -----------
    graphs_data : list
        List of dictionaries containing graph data from generate_undirected_chordal_set
        
    Returns:
    --------
    dict
        Dictionary containing analysis results
    """
    analysis = {}
    
    # Group by edge ratio
    for ratio in set(g['edge_ratio'] for g in graphs_data):
        ratio_graphs = [g for g in graphs_data if g['edge_ratio'] == ratio]
        
        analysis[ratio] = {
            'avg_density': np.mean([g['density'] for g in ratio_graphs]),
            'avg_max_clique_size': np.mean([g['max_clique_size'] for g in ratio_graphs]),
            'avg_num_cliques': np.mean([g['num_cliques'] for g in ratio_graphs]),
            'all_chordal': all(g['is_chordal'] for g in ratio_graphs),
            'all_connected': all(g['is_connected'] for g in ratio_graphs)
        }
    
    return analysis

def visualize_chordal_analysis(graphs_data, output_dir="chordal_analysis"):
    """Create visualizations for chordal graph analysis.
    
    Parameters:
    -----------
    graphs_data : list
        List of dictionaries containing graph data
    output_dir : str
        Directory to save visualization files
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Group by edge ratio
    ratios = sorted(set(g['edge_ratio'] for g in graphs_data))
    
    # Create comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Density distribution
    for ratio in ratios:
        ratio_graphs = [g for g in graphs_data if g['edge_ratio'] == ratio]
        densities = [g['density'] for g in ratio_graphs]
        axes[0, 0].hist(densities, alpha=0.5, label=f'Ratio {ratio:.2f}')
    axes[0, 0].set_xlabel('Graph Density')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Graph Densities')
    axes[0, 0].legend()
    
    # Plot 2: Max clique size distribution
    for ratio in ratios:
        ratio_graphs = [g for g in graphs_data if g['edge_ratio'] == ratio]
        clique_sizes = [g['max_clique_size'] for g in ratio_graphs]
        axes[0, 1].hist(clique_sizes, alpha=0.5, label=f'Ratio {ratio:.2f}')
    axes[0, 1].set_xlabel('Maximum Clique Size')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Distribution of Maximum Clique Sizes')
    axes[0, 1].legend()
    
    # Plot 3: Number of cliques vs edge ratio
    avg_cliques = [np.mean([g['num_cliques'] for g in graphs_data if g['edge_ratio'] == r])
                   for r in ratios]
    axes[1, 0].plot(ratios, avg_cliques, 'o-')
    axes[1, 0].set_xlabel('Edge Ratio')
    axes[1, 0].set_ylabel('Average Number of Cliques')
    axes[1, 0].set_title('Average Number of Cliques vs Edge Ratio')
    
    # Plot 4: Example graph visualization
    example_graph = graphs_data[0]['graph']
    pos = nx.spring_layout(example_graph)
    nx.draw(example_graph, pos, ax=axes[1, 1], with_labels=True,
            node_color='lightblue', node_size=500, font_size=10)
    axes[1, 1].set_title('Example Chordal Graph')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chordal_analysis.png'))
    plt.close()

def genConsistentCbgk(cpdag, num, bgk_type=1):
    """Generate consistent background knowledge.
    
    Args:
        cpdag: A CPDAG as networkx DiGraph
        num: Number of background knowledge constraints to generate
        bgk_type: Type of background knowledge
            1: direct causal relations
            2: non-causal relations
            3: causal relations
            
    Returns:
        List of tuples (x, y, t) where:
            x, y are node labels
            t is the type of relation
    """
    # Get adjacency matrix and node labels
    amat = nx.to_numpy_array(cpdag)
    label_names = list(cpdag.nodes())
    n = len(label_names)
    
    # Convert CPDAG to DAG using pgmpy
    # First get directed and undirected edges
    directed_edges = []
    undirected_edges = []
    for (u, v) in cpdag.edges():
        if (v, u) in cpdag.edges():
            if (u, v) not in undirected_edges and (v, u) not in undirected_edges:
                undirected_edges.append((u, v))
        else:
            directed_edges.append((u, v))
    
    # Create PDAG with these edges
    pdag = PDAG(directed_ebunch=directed_edges, undirected_ebunch=undirected_edges)
    dag = pdag.to_dag()
    dag = nx.DiGraph(dag.edges())  # Convert back to networkx DiGraph
    dag_amat = nx.to_numpy_array(dag)
    
    # Initialize empty background knowledge list
    bgk = []
    
    if bgk_type == 1:  # Direct causal relations
        # Get all directed edges in DAG
        directed_edges = list(dag.edges())
        
        if directed_edges:
            # Randomly sample edges
            num_edges = min(num, len(directed_edges))
            selected_edges = random.sample(directed_edges, num_edges)
            
            # Create background knowledge
            for x, y in selected_edges:
                bgk.append((x, y, 1))
                    
    elif bgk_type == 2:  # Non-causal relations
        # Get all non-adjacent pairs and non-descendant pairs
        non_causal = []
        for i in range(n):
            # Get descendants in DAG
            descendants = nx.descendants(dag, label_names[i])
            non_descendants = set(label_names) - descendants - {label_names[i]}
            
            for j in range(n):
                if i != j:
                    # Add if j is not a descendant of i
                    if label_names[j] in non_descendants:
                        non_causal.append((label_names[i], label_names[j]))
        
        if non_causal:
            # Randomly sample pairs
            num_pairs = min(num, len(non_causal))
            selected_pairs = random.sample(non_causal, num_pairs)
            
            # Create background knowledge
            for x, y in selected_pairs:
                bgk.append((x, y, 0))
                
    elif bgk_type == 3:  # Causal relations
        # Get all directed edges and descendant relations in DAG
        causal_relations = []
        for i in range(n):
            # 使用 networkx 的 descendants 函数找后代
            descendants = nx.descendants(dag, label_names[i])
            for d in descendants:
                causal_relations.append((label_names[i], d))
        
        if causal_relations:
            # Randomly sample relations
            num_relations = min(num, len(causal_relations))
            selected_relations = random.sample(causal_relations, num_relations)
            
            # Create background knowledge
            for x, y in selected_relations:
                bgk.append((x, y, 1))
    
    return bgk

def plot_graph(G, title="", output_path=None, node_labels=None, highlight_edges=None, 
               edge_colors=None, node_colors=None, pos=None, return_pos=False):
    """Create a Graphviz visualization for graphs in causal inference.
    
    Args:
        G: NetworkX graph
        title: Title of the graph
        output_path: Path to save the image
        node_labels: Optional dict mapping nodes to labels
        highlight_edges: Optional set of edges to highlight
        edge_colors: Optional dict mapping edges to colors
        node_colors: Optional dict mapping nodes to colors (新增)
        pos: Optional node positions
        return_pos: Whether to return node positions
    """
    dot = Digraph(comment='Graph')
    dot.attr(rankdir='LR')
    
    # 设置图形属性
    dot.attr('node', shape='circle', style='filled', fillcolor='lightblue',
            width='0.8', height='0.8', fixedsize='true')  # 调整节点大小
    dot.attr('edge', splines='curved')  # 使用曲线边
    
    # 如果没有提供布局，使用更好的布局设置
    if not pos:
        dot.attr(layout='neato',
                overlap='false',     # 避免节点重叠
                sep='+25,25',       # 增加节点间距
                splines='curved',    # 使用曲线
                rankdir='LR')       # 从左到右的布局
    
    # 添加节点
    for node in G.nodes():
        label = str(node_labels[node]) if node_labels and node in node_labels else str(node)
        if pos:
            x, y = pos[node]
            x = (x + 1) * 5
            y = (y + 1) * 5
            # 添加节点颜色支持
            if node_colors and node in node_colors:
                dot.node(str(node), label, pos=f"{x},{y}!", fillcolor=node_colors[node])
            else:
                dot.node(str(node), label, pos=f"{x},{y}!")
        else:
            # 添加节点颜色支持
            if node_colors and node in node_colors:
                dot.node(str(node), label, fillcolor=node_colors[node])
            else:
                dot.node(str(node), label)
    
    # 处理边（保持原有功能不变）
    processed_edges = set()
    for u, v in G.edges():
        if (u, v) not in processed_edges:
            is_undirected = (v, u) in G.edges()
            if is_undirected:
                color = edge_colors.get((u, v)) or edge_colors.get((v, u)) if edge_colors else 'black'
                width = '2.0' if highlight_edges and ((u, v) in highlight_edges or (v, u) in highlight_edges) else '1.5'
                dot.edge(str(u), str(v), 
                        dir='none', 
                        color=color, 
                        penwidth=width)
                processed_edges.add((u, v))
                processed_edges.add((v, u))
            else:
                color = edge_colors.get((u, v)) if edge_colors else 'black'
                width = '2.0' if highlight_edges and (u, v) in highlight_edges else '1.5'
                dot.edge(str(u), str(v), 
                        dir='forward',
                        color=color, 
                        penwidth=width)
                processed_edges.add((u, v))
    
    # 获取布局信息（保持原有功能不变）
    if not pos and return_pos:
        pos = {}
        for line in dot.pipe('dot', format='plain').decode().split('\n'):
            if line.startswith('node'):
                parts = line.split()
                if len(parts) >= 4:
                    node = parts[1]
                    x = float(parts[2]) / 5 - 1  # 转换回[-1,1]范围
                    y = float(parts[3]) / 5 - 1
                    pos[node] = (x, y)
    
    # 生成图像（保持原有功能不变）
    png_data = dot.pipe(format='png')
    encoded_image = base64.b64encode(png_data).decode('utf-8')
    image_data = f'data:image/png;base64,{encoded_image}'
    
    if return_pos:
        return image_data, pos
    return image_data

def plot_dag(G, title="DAG", output_path=None, **kwargs):
    """Plot a Directed Acyclic Graph (DAG).
    
    Args:
        G: NetworkX DiGraph
        title: Title of the graph
        output_path: Path to save the image
        **kwargs: Additional arguments passed to plot_graph
    """
    return plot_graph(G, title=title, output_path=output_path, **kwargs)

def plot_cpdag(G, title="CPDAG", output_path=None, **kwargs):
    """Plot a Completed Partially Directed Acyclic Graph (CPDAG).
    
    Args:
        G: NetworkX DiGraph with some bidirectional edges
        title: Title of the graph
        output_path: Path to save the image
        **kwargs: Additional arguments passed to plot_graph
    """
    return plot_graph(G, title=title, output_path=output_path, **kwargs)

def plot_mpdag(G, title="MPDAG", output_path=None, **kwargs):
    """Plot a Maximally Oriented Partially Directed Acyclic Graph (MPDAG).
    
    Args:
        G: NetworkX DiGraph with some bidirectional edges
        title: Title of the graph
        output_path: Path to save the image
        **kwargs: Additional arguments passed to plot_graph
    """
    return plot_graph(G, title=title, output_path=output_path, **kwargs)

def plot_pdag(G, title="PDAG", output_path=None, **kwargs):
    """Plot a Partially Directed Acyclic Graph (PDAG).
    
    Args:
        G: NetworkX DiGraph with some bidirectional edges
        title: Title of the graph
        output_path: Path to save the image
        **kwargs: Additional arguments passed to plot_graph
    """
    return plot_graph(G, title=title, output_path=output_path, **kwargs)

def plot_graph_with_bgk(G, bgk, title="Graph with Background Knowledge", output_path=None, **kwargs):
    """Plot a graph with background knowledge, showing causal and non-causal edges.
    
    Args:
        G: NetworkX DiGraph
        bgk: List of tuples (x, y, t) where t=1 for causal and t=0 for non-causal
        title: Title of the graph
        output_path: Path to save the image
        **kwargs: Additional arguments passed to plot_graph
    """
    # Create edge colors dictionary for background knowledge
    edge_colors = {}
    for x, y, t in bgk:
        edge_colors[(x, y)] = '#2ecc71' if t == 1 else '#e74c3c'  # Green for causal, red for non-causal
    
    return plot_graph(G, title=title, output_path=output_path, edge_colors=edge_colors, **kwargs)

def plot_graphs_comparison(G1, G2, title1="Graph 1", title2="Graph 2", output_path=None, **kwargs):
    """Plot two graphs side by side for comparison.
    
    Args:
        G1: First NetworkX graph
        G2: Second NetworkX graph
        title1: Title for first graph
        title2: Title for second graph
        output_path: Path to save the image
        **kwargs: Additional arguments passed to plot_graph
    """
    # Create Digraph object
    dot = Digraph()
    
    # Create subgraphs for side-by-side comparison
    with dot.subgraph(name='cluster_0') as c1:
        c1.attr(label=title1)
        # Add nodes and edges for G1
        for node in G1.nodes():
            c1.node(f"1_{node}", str(node))
        
        # Process edges for G1
        processed_edges = set()
        for u, v in G1.edges():
            if (u, v) not in processed_edges:
                is_undirected = (v, u) in G1.edges()
                if is_undirected:
                    c1.edge(f"1_{u}", f"1_{v}", dir='none')
                    processed_edges.add((u, v))
                    processed_edges.add((v, u))
                else:
                    c1.edge(f"1_{u}", f"1_{v}", dir='forward')
                    processed_edges.add((u, v))
    
    with dot.subgraph(name='cluster_1') as c2:
        c2.attr(label=title2)
        # Add nodes and edges for G2
        for node in G2.nodes():
            c2.node(f"2_{node}", str(node))
        
        # Process edges for G2
        processed_edges = set()
        for u, v in G2.edges():
            if (u, v) not in processed_edges:
                is_undirected = (v, u) in G2.edges()
                if is_undirected:
                    c2.edge(f"2_{u}", f"2_{v}", dir='none')
                    processed_edges.add((u, v))
                    processed_edges.add((v, u))
                else:
                    c2.edge(f"2_{u}", f"2_{v}", dir='forward')
                    processed_edges.add((u, v))
    
    # Render and save or return the comparison
    if output_path:
        dot.render(output_path, format='png', cleanup=True)
        return None
    else:
        png_data = dot.pipe(format='png')
        encoded_image = base64.b64encode(png_data).decode('utf-8')
        return f'data:image/png;base64,{encoded_image}'

def plot_weighted_graph(G, title="", output_path=None, edge_colors=None, pos=None, return_pos=False):
    dot = Digraph(comment='Weighted Graph')
    dot.attr(rankdir='LR')
    
    # 设置图形属性
    dot.attr('node', shape='circle', style='filled', fillcolor='lightblue',
            width='0.8', height='0.8', fixedsize='true')  # 调整节点大小
    dot.attr('edge', splines='curved')  # 使用曲线边
    
    # 如果没有提供布局，使用更好的布局设置
    if not pos:
        dot.attr(layout='neato',
                overlap='false',     # 避免节点重叠
                sep='+25,25',       # 增加节点间距
                splines='curved',    # 使用曲线
                rankdir='LR')       # 从左到右的布局
    
    # 添加节点
    for node in G.nodes():
        if pos:
            x, y = pos[node]
            x = (x + 1) * 5
            y = (y + 1) * 5
            dot.node(str(node), str(node), pos=f"{x},{y}!")
        else:
            dot.node(str(node), str(node))
    
    # 处理边
    processed_edges = set()
    for u, v in G.edges():
        if (u, v) not in processed_edges:
            weight = G[u][v].get('weight', 1.0)
            weight_str = f"{abs(weight):.2f}"
            weight_label = f' {weight_str} '
            
            is_undirected = (v, u) in G.edges()
            if is_undirected:
                color = edge_colors.get((u, v)) or edge_colors.get((v, u)) if edge_colors else 'black'
                dot.edge(str(u), str(v), 
                        dir='none', 
                        color=color, 
                        label=weight_label, 
                        fontsize='10',
                        penwidth='1.5')  # 增加边的宽度
                processed_edges.add((u, v))
                processed_edges.add((v, u))
            else:
                color = edge_colors.get((u, v)) if edge_colors else 'black'
                dot.edge(str(u), str(v), 
                        dir='forward', 
                        color=color, 
                        label=weight_label, 
                        fontsize='10',
                        penwidth='1.5')  # 增加边的宽度
                processed_edges.add((u, v))
    
    # 获取布局信息（如果没有提供的话）
    if not pos and return_pos:
        # 使用graphviz的布局信息
        pos = {}
        for line in dot.pipe('dot', format='plain').decode().split('\n'):
            if line.startswith('node'):
                parts = line.split()
                if len(parts) >= 4:
                    node = parts[1]
                    x = float(parts[2]) / 5 - 1  # 转换回[-1,1]范围
                    y = float(parts[3]) / 5 - 1
                    pos[node] = (x, y)
    
    # 生成图像
    png_data = dot.pipe(format='png')
    encoded_image = base64.b64encode(png_data).decode('utf-8')
    image_data = f'data:image/png;base64,{encoded_image}'
    
    if return_pos:
        return image_data, pos
    return image_data