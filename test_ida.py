import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["NUMEXPR_MAX_THREADS"] = "8"
from rcbgk import ida, true_cov
import numpy as np
import networkx as nx

def test_ida():
    # 定义 CPDAG 和 DAG 的邻接矩阵
    amat_cpdag = np.array([[0, 1, 1], 
                          [1, 0, 1], 
                          [1, 1, 0]])  # 无向图
                          
    amat_dag = np.array([[0, 1, 1], 
                         [0, 0, 1], 
                         [0, 0, 0]])   # DAG
    
    # 创建 NetworkX 图对象
    cpdag = nx.DiGraph(amat_cpdag)
    dag = nx.DiGraph(amat_dag)
    
    # 设置测试参数
    x = 0  # 源节点
    y = 2  # 目标节点
    
    # 计算协方差矩阵
    mcov = true_cov(dag)
    
    # 调用 IDA 函数
    effects = ida(x, y, mcov, cpdag)
    
    print(f"CPDAG adjacency matrix:\n{amat_cpdag}")
    print(f"DAG adjacency matrix:\n{amat_dag}")
    print(f"CPDAG: {cpdag}")
    print(f"DAG: {dag}")
    print(f"Covariance matrix:\n{mcov}")
    print(f"Possible causal effects from {x} to {y}: {effects}")

if __name__ == "__main__":
    test_ida()