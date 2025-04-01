library(pcalg)

# 定义 CPDAG 和 DAG 的邻接矩阵
amat.cpdag <- matrix(c(0, 1, 1,
                      1, 0, 1,
                      1, 1, 0), nrow=3, byrow=TRUE)

amat.dag <- matrix(c(0, 1, 1,
                    0, 0, 1,
                    0, 0, 0), nrow=3, byrow=TRUE)

cpdag <- as(amat.cpdag, "graphNEL")
dag <- as(amat.dag, "graphNEL")

# 设置测试参数
x <- 1  # 源节点
y <- 3  # 目标节点

# 计算协方差矩阵
mcov <- trueCov(dag)

# 调用 IDA 函数
effects <- ida(x.pos=x, y.pos=y, mcov=mcov, graphEst=cpdag, method="local", type="cpdag")

# 打印结果
cat("CPDAG adjacency matrix:\n")
print(amat.cpdag)
cat("\nDAG adjacency matrix:\n")
print(amat.dag)
cat("\nCPDAG:\n")
print(cpdag)
cat("\nDAG:\n")
print(dag)
cat("\nCovariance matrix:\n")
print(mcov)
cat("\nPossible causal effects from", x, "to", y, ":", effects, "\n")

