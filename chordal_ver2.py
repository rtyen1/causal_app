from collections import deque
import networkx as nx
import time
##import gmpy2
import random
import math
import pylab as P 
import copy
import pickle
import os 
from sympy import *

global Gback
Gback=[]
global A  
m=symbols('m')

nMax=500
A=[[0 for i in range(nMax)] for j in range(nMax)]    
for j in range(nMax):
    A[0][j]=(-1)**(j) 
    A[j][j]=j+1
    for i in range(1,j):
        A[i][j]=-A[i-1][j]*(j+2-i)/i
        
      

H={0: {1, 2, 3, 4, 5},
 1: {0, 2},
 2: {0, 1},
 3: {0, 4, 5},
 4: {0, 3, 5},
 5: {0, 3, 4}}

Hlist=[{0:{1,2},1:{0,2},2:{0,1}},{3:{4},4:{3}},{5:{}},{6:{}},{7:{8,9,10},8:{7,9,10},9:{7,8,10},10:{7,8,9}}]

H ={0: {1, 2, 3, 4},
 1: {0, 2},
 2: {0, 1,3},
 3: {0, 2,4},
 4: {0, 3}}

H={0: set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), 1: set([0, 2, 5, 6, 7, 8, 9, 10]), 2: set([0, 1, 3, 4, 5, 6, 7, 8, 9, 10]), 3: set([0, 2, 4, 6, 8, 10]), 4: set([0, 2, 3, 5, 6, 8, 9, 10]), 5: set([0, 1, 2, 4, 6, 7, 8, 9, 10]), 6: set([0, 1, 2, 3, 4, 5, 7, 8, 9, 10]), 7: set([0, 1, 2, 5, 6, 8, 9, 10]), 8: set([0, 1, 2, 3, 4, 5, 6, 7, 9, 10]), 9: set([0, 1, 2, 4, 5, 6, 7, 8, 10]), 10: set([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])}
H={0: {1, 2, 3, 4, 5, 6},
 1: {0, 6},
 2: {0, 3, 4, 5},
 3: {0, 2, 4, 5, 6},
 4: {0, 2, 3, 5},
 5: {0, 2, 3, 4, 6},
 6: {0, 1, 3, 5}}
 
H={0: {2, 3, 4, 5, 6},
 1: {4},
 2: {0, 3, 4, 5, 6},
 3: {0, 2, 4, 6},
 4: {0, 1, 2, 3, 5, 6, 7},
 5: {0, 2, 4, 6},
 6: {0, 2, 3, 4, 5},
 7: {4}}


Hlist=[{0:{}},{1:{}}]
Hlist=[{0:set()},{1:{2},2:{1}}]
Hlist123=[{0:{},1:{}},{0:{1},1:{0},2:{}},{0:{},1:{},2:{}},{0:{1},1:{0,2},2:{1,3},3:{2}},{0:{1,2},1:{0,2},2:{0,1},3:{}}]
Hlist4=[{0:{1},1:{0},2:{3},3:{2}},{0:{1,2},1:{0},2:{0},3:{}},{0:{1,2,3},1:{0,3,2},2:{0,1,3},3:{0,1,2},4:{}},{0:{1,2},1:{0,3,2},2:{0,1,3},3:{4,1,2},4:{3}}]


#Htem=map(sepgraph,Hlist123)
#formul1=map(formulsize,Htem)
#Htem=map(sepgraph,Hlist4)
#formul2=map(formulsize,Htem)


#H=random_connected_chordal_ver3(20,130)
#t1=time.time();s1=evalsize(H);tused=time.time()-t1; print(s1,t1)

global dominating
dominating=3
def evalsize(H):    
    ##H  a UCCG in a dictionary with neigbor sets of every vertex
    #The following code obtain the degree of each vertex
    
    pvset=len(H);
    nedge=sum([len(H[i])   for i in H])/2   
    if nedge==pvset-1:     
        return(pvset)         
    if nedge==pvset:
        return(2*pvset) 
    if nedge==pvset*(pvset-1)/2-1:
        return(2*factorial(pvset-1)-factorial(pvset-2))
    if nedge==pvset*(pvset-1)/2:
        return(factorial(pvset))
    if nedge==pvset*(pvset-1)/2-2:
        return((pvset**2-pvset-4)*factorial(pvset-3))
  
 
    fulld=set([i  for i in H if len(H[i])==pvset-1])
    nf=len(fulld)                   
    ##Below, reduce the base graph H if possible 
    if nf>dominating:
        #        print "formul"
        Hlist={}
        for index in H:
            if index not in fulld:
                Hlist[index]=  H[index]-fulld 
        Hlist=sepgraph(Hlist)    # need to be revised
        formul=formulsize(Hlist)
        return(formul.subs(m,nf))
    else:
        num1=0
        for kkk in H:
            post=postrootv3(kkk,H)
            #print(post)
            temcoun=1
            for temneighbor in post:       
                temcoun*=evalsize(temneighbor)
            num1+=temcoun
        #print(6)
##        print("degree not 1:", num1)
        return(num1)

def evalHlist(Hlist):
    size=1
    for H in Hlist:
        size=size*evalsize(H)
    return size
        

def formulsize(Hlist):
    ##deduce the formula of the size of the graph obtained by extended H by m full degree vertices
    ## H has J UCCG, each in a dictrionary of a list, looks like [{},{}]
    pn=[[len(i),sum([len(i[j]) for j in i])/2] for i in Hlist]   
    J=len(pn)      
    #empty graph
    if J==0:
        return(factorial(m))
    if J==1:
        ##Tree
        if pn[0][0]==pn[0][1]+1:            
            return((pn[0][1]*m**2+(2*pn[0][0]-1)*m+pn[0][0])*factorial(m))
        ## Tree plus 1
        if pn[0][0]==pn[0][1]:
            return(( m**3+2*pn[0][0]*m**2 +(4*pn[0][0]-1)*m+2*pn[0][0])*factorial(m))        

       ## full degree vertices
        fulld=set([i  for i in Hlist[0] if len(Hlist[0][i])==len(Hlist[0])-1])
        nf=len(fulld)  
        if nf>0:
            Hlisttem={}     
            for index in Hlist[0]:
                if index not in fulld:
                    Hlisttem[index] =  Hlist[0][index]-fulld 
            Hlisttem=sepgraph(Hlisttem)  
            formutem=formulsize(Hlisttem)
            #formutem=expand()
            #formulbase=[Hlist[0],formutem]
            #Gback.append(formulbase)
            return formutem.replace(m,m+nf)
        
    ## isolated edges
    if all([(pni[0]==2) & (pni[1]==1) for pni in pn]):     
        return(2**(J-1)*(J*m**2+3*J*m+2)*factorial(m))
    
    ## isolated vertex
    non_isovertex=[len(i) for i in Hlist if len(i)>1]     
    k=len(non_isovertex)
    if k<J:
        non_isovertexH=[i for i in Hlist if len(i)>1]
        res=formulsize(non_isovertexH)+evalHlist(non_isovertexH)*(J-k)*m*factorial(m)
        #res=expand(res)
        return(res)

 
            
 
    ###cursive
    gm=increment(Hlist) 
    alpha=gm.all_coeffs()
    alpha.reverse()
    d=len(alpha)     
    ##d is degree of f; d-1 is degree of g
    beta=[0 for i in range(d)]    
    beta[d-1]=alpha[d-1]/A[d-1][d-1]
    for i in range(d-2,-1,-1):
        ss=0        
        for k in range(i+1,d):
            ss=ss+beta[k]*A[i][k]
        beta[i]=(alpha[i]-ss)/A[i][i]   
    beta.reverse() 
    beta.append(evalHlist(Hlist)) 
    res=Poly(beta,m)*factorial(m)
    #formulbase=[Hlist,res]
    #Gback.append(formulbase)
    return(res)

def coregraph(H):
    de =set()
    p=len(H)
    for i in H:
        if len(H[i])<p-1:
            de.add(i)
            
    coreg={}
    for i in H:
        if i in de:
            coreg[i]=H[i].intersection(de)
    return(coreg)
  
def extgraph(H,k):
    G={}
    if isinstance(H,list): 
        for items in H:
            G.update(items)
    else:
        G.update(H)
    Key=set(G.keys() )
    Maxx=max(Key)    
    M=set(range(Maxx+1,Maxx+k+1))
    Allkey=M.union(Key)
    for items in G:
        G[items] 
        G[items].update(M)
    for k in M:
        G[k]=Allkey.difference(set([k]))
    return(G)
    
def nbplotback(H,pos=None,with_lab=True):
    G=[]
    if isinstance(H,list): 
        for items in H:
            for nodei in items:
                for nodej in items[nodei]:
                    if nodej>nodei:
                        G.append([nodei,nodej])
           
    else:
        for nodei in H:
            for nodej in H[nodei]:
                if nodej>nodei:
                    G.append([nodei,nodej])
                      
    nxg=nx.Graph(G) 
    
    if pos is None:
        pos = nx.spring_layout(nxg)    
    nx.draw(nxg,pos=pos,with_labels=with_lab,node_color="w") 
    return(pos)
            
                        
def nbplot(H,pos=None,with_lab=True,mywidth=2):    
    if isinstance(H,list): 
        G={}
        for items in H:
            G.update(items) 
    else:
        G=H
                      
    nxg=nx.Graph(G) 
    
    if pos is None:
        pos = nx.spring_layout(nxg)    
    nx.draw(nxg,pos=pos,with_labels=with_lab,node_color="w",width=mywidth) 
    return(pos)
           
    
def sepgraph(G):
    left=set(G.keys()) 
    H=[]
    while len(left)>0:
        current=left.pop()
        visited={current}
        tovisit=G[current].copy()
        Tg={}
        Tg[current]=G[current].copy()        
        while len(tovisit)>0:
            tem=tovisit.pop() 
            visited.add(tem)
            Tg[tem]=G[tem].copy()
            tovisit.update(G[tem].difference(visited))
            left.remove(tem)
        H.append(Tg)

    return(H)
            
            
  
        
def increment(Hlist):  
    #2016-4-30 use combsimp for simplify factorial(m+100)/factorial(m)
    J=len(Hlist)
    s=0;c=1
    ss=[];cc=[]
    for j in range(J):
        s1=0;c1=0
        for v in Hlist[j]:
            
            ##Hv: chain components of the v-rooted graph of Hj, 
            Hv=postrootv3(v,Hlist[j])
            ## HNv the subgraph of Hj over Nv, including the isolated vertices
            HNv=[]
            
            ## the neighbor set of v in Hj
            Nv={nv for nv in Hlist[j][v]}    
            
            for subt in Hv:
                if subt.keys()[0] in Nv:
                    HNv.append(subt)
                    Nv=Nv- set(subt.keys())                
            for vv in Nv:
                HNv.append({vv:set([])})
            
            sjnv=evalHlist(HNv)
            sjv=evalHlist(Hv)
            s1=s1+formulsize(HNv)*sjv/sjnv
            c1=c1+sjv
        ss.append(s1)
        cc.append(c1)
        c=c*c1
    for kkk in range(len(ss)):        
        s=s+ss[kkk]*c/cc[kkk] 
    s=combsimp(s/factorial(m))  
    s=Poly(s, m)   
    return(s)
        
                
  
    


def random_connected_chordal_ver1(nnodes,nedges):
    ##generate an undirected and  connected chordal graph
    ##use  package networkx to check chordal graph
    ##inefficent
    ##error for nedges is large
    
    ung={i:set([]) for i in range(nnodes)}
    left=set(range(1,nnodes))
    connected=[0]
    for i in range(nnodes-1):
        x=random.sample(connected,1)[0]
        y=random.sample(left,1)[0]
        ung[x].add(y)
        ung[y].add(x)
        connected.append(y)
        left.remove(y)

    leftedn=nedges-nnodes+1
    nodes=range(nnodes)
    
    for i in range(leftedn):
        test=1
        while test:
            test1=1
            while test1:
                u=random.choice(nodes)
                ulist=list(ung[u])
                le=len(ulist)
                ## all pairs that are not adjancent in neigbors of u 
                unedgebase=[[ulist[ii],ulist[j]] for ii in range(le) for j in range(ii+1,le) if not(ulist[ii] in ung[ulist[j]])]
                if len(unedgebase)>0:
                    test1=False
            ed=random.choice(unedgebase)
            
            ung[ed[0]].add(ed[1])
            ung[ed[1]].add(ed[0])
            e=[(iii,j) for iii in ung for j in ung[iii] if i<j]
            G=nx.Graph(e)
            if nx.is_chordal(G):        
                test=False
            else:
                ung[ed[0]].remove(ed[1])
                ung[ed[1]].remove(ed[0])
    return(ung)

def sub_graph(G,subset):
    return {i:G[i]&subset for i in subset}
def max_deg(G):
    temp=0
    ind=-1
    for i in G:
        if len(G[i])>temp:
            temp=len(G[i])
            ind=i
    return ind

def is_completed(G):
    nn=len(G)
    nedge=sum([len(G[i]) for i in G])/2
    if nedge==nn*(nn-1)/2:
        return True
    else:
        return False
def perfect_order(G):
    ##relabel vertices such that 1,2,...p is a perfect order
    nodes=G.keys()
    order=[random.choice(nodes)]
    nodes.remove(order[0])
    while len(nodes)>0:
        tem1=sub_graph(G,set(nodes))         
        rem1=max_deg(tem1)
        if rem1==-1:
            order.extend(nodes)
            nodes=[]
        else:
            nodes.remove(rem1)
            order.append(rem1)
    order_index={i:order.index(i) for i in order}
    orderung={order_index[i]:set([order_index[j] for j in G[i]]) for i in G}
    return orderung

def perfect_non_adj(G):
    temp=set()
    for i in G:
        ulist=list(G[i])
        le=len(ulist)
        for i in range(le):
            for j in range(i+1,le):
                if not(ulist[i] in G[ulist[j]]):
                    mj=max(ulist[i],ulist[j])
                    mi=min(ulist[i],ulist[j])
                    subsetmj=set([mjj for mjj in G[mj] if mjj <mj])
                    if subsetmj.issubset(G[mi]):
                        temp.update([(ulist[i],ulist[j])])    
    return temp

def perfect_non_adj_ver2(G):
    #list all egdes
    temp=set()
    nodes=list(G.keys())  # 修复 dict.keys()
    random.shuffle(nodes)
    for i in nodes:
        ulist=list(G[i])  # 修复 dict.values()
        random.shuffle(ulist)
        le=len(ulist)
        for i in range(le):
            for j in range(i+1,le):
                if not(ulist[i] in G[ulist[j]]):
                    G[ulist[i]].add(ulist[j])
                    G[ulist[j]].add(ulist[i])
                    e=[(ii,jj) for ii in G for jj in G[ii] if ii<jj]
                    G1=nx.Graph(e)            
                    if nx.is_chordal(G1):                        
                        temp.update([(ulist[i],ulist[j])])
                    G[ulist[i]].remove(ulist[j])
                    G[ulist[j]].remove(ulist[i])
    return temp 
    
def perfect_non_adj_ver3(G):
    #add an edge directedly, randomly
    #output a graph 
    nodes=list(G.keys())  # 修复 dict.keys()
    random.shuffle(nodes)
    for i in nodes:
        ulist=list(G[i])  # 修复 dict.values()
        random.shuffle(ulist)
        le=len(ulist)
        for i in range(le):
            for j in range(i+1,le):
                if not(ulist[i] in G[ulist[j]]):
                    G[ulist[i]].add(ulist[j])
                    G[ulist[j]].add(ulist[i])
                    e=[(ii,jj) for ii in G for jj in G[ii] if ii<jj]
                    G1=nx.Graph(e)            
                    if nx.is_chordal(G1):                        
                        return G
                    G[ulist[i]].remove(ulist[j])
                    G[ulist[j]].remove(ulist[i])
    return -1

def addEdge(G):
    #add an edge directedly, randomly
    #output a graph 
    e=[(ii,jj) for ii in G for jj in G[ii] if ii<jj]
    nodes=list(G.keys())  # 修复 dict.keys()
    random.shuffle(nodes)
    for i in nodes:
        jset=[j for j in nodes if j>i]
        random.shuffle(jset)
        for j in jset:
            if not (j in G[i]):                
                e.append((i,j))
                G1=nx.Graph(e)
                if nx.is_chordal(G1):
                    G[i].add(j)
                    G[j].add(i)
                    return G
                e.pop()                
    return -1
    
    
def tree_construct(p):   
    
    ung={i:set([]) for i in range(p)}
    left=set(range(1,p))
    connected=[0]
    for i in range(p-1):
        x=random.sample(connected,1)[0]
        y=random.sample(left,1)[0]
        ung[x].add(y)
        ung[y].add(x)
        connected.append(y)
        left.remove(y)
    return ung

def gen_conn_chordal(nnodes,nedges):   ## 
    ung=tree_construct(nnodes)
    leftedn=nedges-nnodes+1
    for i in range(leftedn):
        ung=addEdge(ung)
    return(ung)

def random_connected_chordal(nnodes,nedges):
    ##generate an undirected and  connected chordal graph 
    ## count the size 
    ##randomly choose one from above 
    ##2015-4-18
    
    ung=tree_construct(nnodes)
    #above generate a tree
    res=[]
    leftedn=nedges-nnodes+1
    for i in range(leftedn):
        ung=addEdge(ung)
        print(i)
        if i> 14*nnodes:
            t1=time.time();x2=count(ung,i);tused2=time.time()-t1; 
            t1=time.time();x1=evalsize(ung);tused1=time.time()-t1; 
            res.append([i,[x2,tused2],[x1,tused1]])
            print(i,[x2,tused2,tused1,tused2-tused1])
            
        if ung==-1:
            print("error")
            return ung
    return(res)  
    
def random_chordal(nnodes,nedges):
    ##generate an undirected and  connected chordal graph 
    ##randomly choose one from above 
    ##2015-4-18
    res=[]
    ung={i:set() for i in range(nnodes)} 
    leftedn=nedges
    for i in range(leftedn):
        ung=addEdge(ung)
        t1=time.time();x2=evalsize(ung);tused2=time.time()-t1; 
        res.append([i,x2,tused2])
        print(i)
        if ung==-1:
            print("error")
            return ung
    return(ung)  

def random_connected_chordal_ver3(nnodes,nedges):
    ##generate an undirected and  connected chordal graph
 
    ##randomly choose one from above 
    
    ung={i:set([]) for i in range(nnodes)}
    left=set(range(1,nnodes))
    connected=[0]
    for i in range(nnodes-1):
        x=random.sample(connected,1)[0]
        y=random.sample(list(left),1)[0]
        ung[x].add(y)
        ung[y].add(x)
        connected.append(y)
        left.remove(y)
    #above generate a tree
 
    leftedn=nedges-nnodes+1
    for i in range(leftedn):
        ung=perfect_non_adj_ver3(ung)
        if ung==-1:
            return ung
    return(ung)    



def random_connected_chordal_ver2(nnodes,nedges):
    ##generate an undirected and  connected chordal graph
 
    ##randomly choose one from above 
    
    ung={i:set([]) for i in range(nnodes)}
    left=set(range(1,nnodes))
    connected=[0]
    for i in range(nnodes-1):
        x=random.sample(connected,1)[0]
        y=random.sample(left,1)[0]
        ung[x].add(y)
        ung[y].add(x)
        connected.append(y)
        left.remove(y)
    #above generate a tree
    ung=perfect_order(ung)
    ##above generate order tree,1,2,p is a perfect order.

    leftedn=nedges-nnodes+1
    for i in range(leftedn):
        nadj=perfect_non_adj_ver2(ung)
        if len(nadj)==0:
            print("Error")
            return ung
        ##above generate all edge can add        
        addedge=random.sample(nadj,1)[0]
        ung[addedge[0]].add(addedge[1])
        ung[addedge[1]].add(addedge[0])
        ung=perfect_order(ung)
   
    return(ung)
##ung=ranconnchordal(50,50)      
##e=[(i,j) for i in ung for j in ung[i] if i<j]
##G=nx.Graph(e)            
##nx.is_chordal(G)


##        
##        
##        
##        
##     
##test1pdag8={}
##test1pdag8[1]=set([ 5])
##test1pdag8[2]=set([ 3,4,5])
##test1pdag8[3]=set([2, 4,5])
##test1pdag8[4]=set([ 2,3,5])
##test1pdag8[5]=set([1,2,3,4])
##


def ranconngrap(nnodes,nedges):
    #this function generate a connected undirected chordal graph
    ## first generate a tree
    ## randomly choose an edges and  a adjacent edge of it, construct a trangle 
 
    left=set(range(1,nnodes))
    connected=[0]
    dag=[[] for i in range(nedges)]
    for i in range(nnodes-1):
        x=random.sample(connected,1)[0]
        y=random.sample(left,1)[0]
        connected.append(y)
        left.remove(y)
        if x<y:
            dag[i]=[x,y]
        else:
            dag[i]=[y,x]
    j=nnodes-1
    vset=range(nnodes)
    while j<nedges:
        tem=random.sample(vset,2)
        if tem[0]>tem[1]:
            tem.reverse()
        if tem not in dag:
            dag[j]=tem[:]
            j+=1;
    return(dag)
            
            
        
def rangrap(nnodes,nedges):
    #dag without constraints
    nodeorder=random.sample(range(nnodes),nnodes)
    index=[[i,j] for i in range(nnodes) for j in range (i+1,nnodes)]
    edgesrandom=random.sample(index,nedges)
    for i in range(nedges):
        edgesrandom[i]=[nodeorder[edgesrandom[i][0]],nodeorder[edgesrandom[i][1]]]
    return(edgesrandom)
                
        

def tporder(dag,p):
    #algorithm from BNT
    desdag=[[] for i in range(p)]
    parnum=[0 for i in range(p)]
    NEdge=len(dag)
    for i in range(NEdge):
        desdag[dag[i][0]].append(dag[i][1])
        parnum[dag[i][1]]+=1
    zero_indeg=deque([i for i in range(p) if parnum[i]==0])
    s=0
    tpord=[0 for i in range(p)]
    while len(zero_indeg)>0:
        v=zero_indeg.popleft()
        tpord[s]=v
        s+=1
        cs=desdag[v]
        if len(cs)==0:
            ##next
            continue
        for j in cs:
            parnum[j]-=1
            if parnum[j]==0:
                zero_indeg.appendleft(j)
    return(tpord)
    
def edgeorder(dag,p):
    #algorithm from chickering
##    dag1=copy.deepcopy(dag)
    n=len(dag)
    tpord=tporder(dag,p)
    nordmap=[0 for i in range(p)]
    for i in range(p):
        nordmap[tpord[i]]=i        

    a=sorted(dag,key=lambda x:(nordmap[x[1]], -nordmap[x[0]]))
   
    return(a)
              
def labeledge(dag1,p):
    #
    #compelled 1; reversible=-1; unkown=0
##    dag1=dag[:]
    ansdag=[[] for i in range(p)]
    edgenum=[[] for i in range(p)]
    edgeord=edgeorder(dag1,p)
    NEdge=len(dag1)
    for i in range(NEdge):
        ansdag[edgeord[i][1]].append(edgeord[i][0])
        edgenum[edgeord[i][1]].append(i)
    unlabled=len(dag1)
    edgelabel=[0 for i in range(unlabled)]
    while unlabled>0:
        v=edgelabel.index(0) #step 4 , v is index of x to y
        x=edgeord[v][0]
        y=edgeord[v][1]
        for wi in [i for i in range(len(ansdag[x])) if edgelabel[edgenum[x][i]]==1]:
            if ansdag[x][wi] not in ansdag[y]:
                for j in edgenum[y]:
                    if edgelabel[j]==0:
                        unlabled-=1
                    edgelabel[j]=1
                continue
            else:
                wyi=ansdag[y].index(ansdag[x][wi])
                if edgelabel[edgenum[y][wyi]]==0:
                    unlabled-=1
                edgelabel[edgenum[y][wyi]]=1
        #begin step 8         
        if len([z for z in ansdag[y] if z!=x and z not in ansdag[x]])>0:
            for yi in edgenum[y]:
                if edgelabel[yi]==0:
                    unlabled-=1
                    edgelabel[yi]=1
        else:
            for yi in edgenum[y]:
                if edgelabel[yi]==0:
                    unlabled-=1
                    edgelabel[yi]=-1
    for i in range(len(dag1)):
        edgeord[i].append(edgelabel[i])            
    return(edgeord)

def cceg(dag,p):
    # generate chain components in the essential graph
    labg=labeledge(dag,p)
    chaincomps=[]
    diedges=[]
    for i in range(len(dag)):
        if labg[i][2]==-1:
            chaincomps.append([labg[i][0],labg[i][1]])
        else:
            diedges.append([labg[i][0],labg[i][1]])
##    out=conngrap(chaincomps)
##    return(out)
    neighbor={}
    #neighbors in eg graph, all nodes in chain component
    fatherset={}
    #father set in eg graph, some children are in chain component
    for i in range(len(chaincomps)):
        if not (neighbor.has_key(chaincomps[i][1])):
            neighbor[chaincomps[i][1]]=set([chaincomps[i][0]])
        else:
            neighbor[chaincomps[i][1]].add(chaincomps[i][0])
        if not (neighbor.has_key(chaincomps[i][0])):
            neighbor[chaincomps[i][0]]=set([chaincomps[i][1]])
        else:
            neighbor[chaincomps[i][0]].add(chaincomps[i][1])
    for i in range(len(diedges)):
        if not (fatherset.has_key(diedges[i][1])):
            fatherset[diedges[i][1]]=set([diedges[i][0]])
        else:
            fatherset[diedges[i][1]].add(diedges[i][0])
    leftnodes=set(neighbor.keys())
    ccc=[]
    while len(leftnodes)>0:
        visited=set()
        to_visit=[leftnodes.pop()]
        while len(to_visit)!= 0:
            v=to_visit.pop()
            visited.add(v)
            if v in leftnodes:
                leftnodes.remove(v)
            to_visit.extend(neighbor[v]&leftnodes)
        ccc.append(visited.copy())
    ccn=len(ccc)
    ccedgen=[0 for i in range(ccn)]
    neighborset=[{} for i in range(ccn)]
    node2cc={}
    for i in range(ccn):
        for nodesss in ccc[i]:
            neighborset[i][nodesss]=neighbor[nodesss]
            node2cc[nodesss]=i
            ccedgen[i]+=len(neighbor[nodesss])
    ccedgen=map(lambda a: a/2,ccedgen)
    output={}
    output["Neiset"]=neighborset
    
    output["nodes2neiset"]=node2cc
    output["diedges"]=diedges
    output["fatherset"]=fatherset
    output["ccedgeno"]=ccedgen
    output["undiedges"]=chaincomps
    
    return(output)


#def findfather(vset,pdag):

def egdataextend(pdag):
    neighborset={}
    for i in pdag['Neiset']:
        for j in list(i.keys()):  # 修复 dict.keys()
            neighborset[j]=copy.copy(i[j])
    sonset=findsonset(pdag)
    adjset=findadj(pdag)
    fatherandnei=find_fandn(pdag)
    pdag['neighborset']=neighborset
    pdag['sonset']=sonset
    pdag['adjset']=adjset
    pdag['fatherandnei']=fatherandnei
    return pdag
    


def findneighbor1(vset,neigh):
    #only used in count function
    neighbor={}
    for i in vset:
        neighbor[i]=neigh[i]&vset
    return(neighbor)

    
        
## dictionary is a slow data structure for index.
## another test show dictionary is faster?

def postrootv2(root,Neighbor):
    ##vset is set type data
    ##algorithm to get a constrainted essential graph with root as out-point node.
    ##in counting functoin calculate neighbor out of postroot function
    vset=set(Neighbor.keys())  # 修复 dict.keys()
    if root not in vset:
        print("root is not a node of subgraph")            
    currentnode=set([root])
    ccc=list()
    left=vset-currentnode
    while len(left)!=0:
        res=set()
        for i in currentnode:
            res|=Neighbor[i]
        res&=left.copy()
        resback=res.copy()
        resneighbor={}
        while len(res)>0:
            aa=res.pop()            
            for jj in (Neighbor[aa] & res):
                if (Neighbor[aa]&currentnode)==(Neighbor[jj]&currentnode):
                    if aa in resneighbor:  # 修复 has_key
                        resneighbor[aa].add(jj)
                    else:
                        resneighbor[aa]=set([jj])
                    if jj in resneighbor:  # 修复 has_key
                        resneighbor[jj].add(aa)
                    else:
                        resneighbor[jj]=set([aa])
                    
 
 
        ccc.extend(sepgraph(resneighbor))

        
        currentnode=resback.copy()
        left-=resback.copy()
        


    return(ccc)
    

def postrootv3(root,Neighbor):
    ##vset is set type data
    ##algorithm to get a constrainted essential graph with root as out-point node.
    ##in counting functoin calculate neighbor out of postroot function
    vset=set(Neighbor.keys())  # 修复 dict.keys()
    if root not in vset:
        print("root is not a node of subgraph")            
    currentnode=set([root])
    ccc=list()
    left=vset-currentnode
    while len(left)!=0:
        res=set()
        fatherdict={}
        for i in currentnode:
            nodestem=Neighbor[i].intersection(left)            
            res|= nodestem
            for nodesi in nodestem:                
                if nodesi in fatherdict:  # 修复 has_key
                    fatherdict[nodesi].add(i)
                else:
                    fatherdict[nodesi]=set([i])            
        resg={i:Neighbor[i]&res for i in res}
        Do=True
        while Do:
            Do=False            
            for i in list(resg.keys()):  # 修复 dict.keys()
                removetem=[]
                for j in resg[i]:
                    if not fatherdict[i].issubset(Neighbor[j]):
                        fatherdict[j].add(i)
                        removetem.append((i,j))
                        Do=True
                for edge in removetem:
                    resg[edge[0]].remove(edge[1])
                    resg[edge[1]].remove(edge[0])
            resg={i:resg[i] for i in resg if len(resg[i])>0}
        ccc.extend(sepgraph(resg))   
        currentnode=res
        left-=res
    return(ccc)
        
        
def postroot(root,Neighbor):
    ##vset is set type data
    ##algorithm to get a constrainted essential graph with root as out-point node.
    ##in counting functoin calculate neighbor out of postroot function
    vset=set(Neighbor.keys())
    if root not in vset:
        print("root is not a node of subgraph")            
    currentnode=set([root])
    ccc=list()
    left=vset-currentnode
    while len(left)!=0:
        res=set()
        fatherdict={}
        for i in currentnode:
            nodestem=Neighbor[i].intersection(left)            
            res|= nodestem
            for nodesi in nodestem:                
                if fatherdict.has_key(nodesi):
                    fatherdict[nodesi].add(i)
                else:
                    fatherdict[nodesi]=set([i])            
        resg={i:Neighbor[i]&res for i in res}
        Do=True
        while Do:
            Do=False            
            for i in  resg:      
                removetem=[]
                for j in resg[i]:
                    if not fatherdict[i].issubset(Neighbor[j]):
                        fatherdict[j].add(i)
                        removetem.append((i,j))
                        Do=True
                for edge in removetem:
                    resg[edge[0]].remove(edge[1])
                    resg[edge[1]].remove(edge[0])
            resg={i:resg[i] for i in resg if len(resg[i])>0 }
        ccc.extend(sepgraph(resg))   
        currentnode=res
        left-=res
        
    ccc=[item.keys() for item in ccc]
    nccc=len(ccc)
    ccnodesn=[0 for i in range(nccc)]
    ccedgesn=[0 for i in range(nccc)]
    for i in range(nccc):
        ccnodesn[i]=len(ccc[i])
        tem=set([j for j in ccc[i]])
        while len(tem)!=0:
            ss=tem.pop()
            ccedgesn[i]+=len(Neighbor[ss]&tem)
    return(ccc,ccnodesn,ccedgesn)

#aa,bb,cc = postroot(root,vset)
#d=time.time()-a
#print(d)
#def countec(cc):


def count(neighbor,nedge):    
    #recursive version of count markov equivlence class, there could be too many
    #nested recursive steps to run normally
    nedge=sum([len(neighbor[i])   for i in neighbor])/2  
    p=len(neighbor)
    vset=neighbor.keys()
    ##        print("case 1:",p)
    if nedge==p-1: 
        return(p)
    
    if nedge==p:
        #print(1)
       ##    print("case 2:",2*p)        

        return(2*p)   

    if nedge==p*(p-1)/2-1:
        #print(2)
       ##        print("case 3:",2*myfac(p-1)-myfac(p-2))
        return(2*factorial(p-1)-factorial(p-2))
    if nedge==p*(p-1)/2:
        #print(vset)
        #print(3)
        ##        print("case 4:",myfac(p))
        return(factorial(p))
    if nedge==p*(p-1)/2-2:
        return((p*p-p-4)*factorial(p-3))
    #    k=0
    #    for i in vset:
    #        if k==0:
    #            minnode=i
    #            tem=len(neighbor[i])
    #            mindeg=tem
    #            k+=1
    #        elif tem>len(neighbor[i]):
    #            minnode=i
    #            mindeg=len(neighbor[i])
    #            tem=mindeg
    #            k+=1      
    #    if mindeg==1:
    ###        print("degree=1")
    #
    #  
    #        post,cnp,cne=postroot(minnode,neighbor)
    #        #print(post)
    #        temcoun1=1
    #        for kk in range(len(post)):
    #            temneighbor=findneighbor1(post[kk],neighbor)
    #            temcoun1*=count(temneighbor,cne[kk])
    #        temnei=neighbor.copy()
    #        del temnei[minnode]
    #        
    #        temcoun2=count(temnei,nedge-1)
    #        #print(5)
    ###        print("degree 1:",temcoun1+temcoun2)
    #        return(temcoun1+temcoun2)
    #    else:
            #other cases with each nodes as root.
     
    ##        print("degree!=1")
    num1=0
    for kkk in vset:
        post,cnp,cne=postroot(kkk,neighbor)
        #print(post)
        temcoun=1
        for kk in range(len(post)):
            temneighbor={i:neighbor[i]&set(post[kk]) for i in post[kk]}
            temcoun*=count(temneighbor,cne[kk])
        num1+=temcoun
    #print(6)
     ##      print("degree not 1:", num1)
    return(num1)
    
##    



global  gcall
global  gdeepth

def count1(neighbor,nedge):
    #recursive version of count markov equivlence class, there could be too many
    #nested recursive steps to run normally
    #to record the times of calling funtion count1 and the depth of recursion,
    #but this depth might not be the real depth
    global gcall
    global gdeepth
    gcall+=1
##    print("call me")

    p=len(neighbor)
    vset=neighbor.keys()
    if nedge==p-1:
        #print(vset)
        #print(0)
##        print("case 1:",p)        
        return(p)
    
    
    if nedge==p:
        #print(1)
##        print("case 2:",2*p)        

        return(2*p)
    

    if nedge==p*(p-1)/2-1:
        #print(2)
##        print("case 3:",2*myfac(p-1)-myfac(p-2))
        return(2*factorial(p-1)-factorial(p-2))
    if nedge==p*(p-1)/2:
        #print(vset)
        #print(3)
##        print("case 4:",myfac(p))
        return(factorial(p))
    if nedge==p*(p-1)/2-2:
        return((p*p-p-4)*factorial(p-3))
    k=0
    for i in vset:
        if k==0:
            minnode=i
            tem=len(neighbor[i])
            mindeg=tem
            k+=1
        elif tem>len(neighbor[i]):
            minnode=i
            mindeg=len(neighbor[i])
            tem=mindeg
            k+=1      
    if mindeg==1:
##        print("degree=1")

        gdeepth+=1
        post,cnp,cne=postroot(minnode,neighbor)
        #print(post)
        temcoun1=1
        for kk in range(len(post)):
            temneighbor={i:neighbor[i]&set(post[kk]) for i in post[kk]}
            temcoun1*=count1(temneighbor,cne[kk])
        temnei=neighbor.copy()
        del temnei[minnode]
        
        temcoun2=count(temnei,nedge-1)
        #print(5)
##        print("degree 1:",temcoun1+temcoun2)
        return(temcoun1+temcoun2)
    else:
        #other cases with each nodes as root.
        gdeepth+=1
##        print("degree!=1")
        num1=0
        for kkk in vset:
            post,cnp,cne=postroot(kkk,neighbor)
            #print(post)
            temcoun=1
            for kk in range(len(post)):
                temneighbor={i:neighbor[i]&set(post[kk]) for i in post[kk]}
                temcoun*=count1(temneighbor,cne[kk])
            num1+=temcoun
        #print(6)
##        print("degree not 1:", num1)
        return(num1)
    
##    
##


global  gcall
global  gdepth

gcall=0
gdepth=0
def count2(neighbor,nedge,dep):
    #recursive version of count markov equivlence class, there could be too many
    #nested recursive steps to run normally
    #to record the times of calling funtion count1 and the depth of recursion,
    #in this version, i modify the method to count depth of recursion.
    # an addition argment dep to record the depth, we pass the max to gdepth
    global gcall
    global gdepth
    gcall+=1
    if gdepth<dep:
        gdepth=dep
##        print dep
##    print("call me")

    p=len(neighbor)
    vset=neighbor.keys()
    if nedge==p-1:
        return(p)    
    if nedge==p:        
        return(2*p)
    if nedge==p*(p-1)/2-1:
        return(2*myfac(p-1)-myfac(p-2))
    if nedge==p*(p-1)/2:
        return(myfac(p))
    if nedge==p*(p-1)/2-2:
        return((p*p-p-4)*myfac(p-3))
    k=0
    for i in vset:
        if k==0:
            minnode=i
            tem=len(neighbor[i])
            mindeg=tem
            k+=1
        elif tem>len(neighbor[i]):
            minnode=i
            mindeg=len(neighbor[i])
            tem=mindeg
            k+=1      
    if mindeg==1:
        post,cnp,cne=postroot(minnode,neighbor)
                
        #print(post)
        temcoun1=1
        for kk in range(len(post)):
            temneighbor={i:neighbor[i]&set(post[kk]) for i in post[kk]}
            temcoun1*=count2(temneighbor,cne[kk],dep+1)
        temnei=neighbor.copy()
        del temnei[minnode]
        
        temcoun2=count2(temnei,nedge-1,dep)
        return(temcoun1+temcoun2)
    else:        
        num1=0
        for kkk in vset:
            post,cnp,cne=postroot(kkk,neighbor)            
            temcoun=1
            for kk in range(len(post)):
                temneighbor={i:neighbor[i]&set(post[kk]) for i in post[kk]}
                temcoun*=count2(temneighbor,cne[kk],dep+1)
            num1+=temcoun
        return(num1)


def countpdag(pdag):
    fatherset=pdag['fatherset']
    neighborsets=pdag['Neiset']
    noedges=pdag['ccedgeno']
    out=1
    cnn=len(neighborsets)
    for ii in range(cnn):
        out*=count(neighborsets[ii],int(noedges[ii]))
    return out


##
        
def countsci(neighbor,nedge):
    num=count(neighbor,nedge)

def myfac(n):
    k=1
    for i in range(1,n+1):
        k*=i
    return(k)


##
##def simulation1(p,ne):    
##    dag=rangrap(p,ne)
##    pdag=cceg(dag,p)
##    n=countpdag(pdag)
##    print(n)

##
##p=50
##gcalltemp=[]
##gdeepthtemp=[]
##for i in range(1):

##@profile
def mysimu(p,n):
    global gcall
    global gdeepth
    t1=time.time()
    ung=random_connected_chordal_ver3(p,n)
    t1=time.time()-t1
##    print time.time()-t1
    if not(ung==-1):
        e=[(i,j) for i in ung for j in ung[i] if i<j]
        G=nx.Graph(e)
##        print nx.is_chordal(G)

        t2=time.time()
        a=count(ung,n)
        t2=time.time()-t2
##        print time.time()-t2
##        print a
##        nx.draw(G)
##        P.draw()
##        P.show()
    #size, time to create chordal graph, time to count, times to call
    return [a,t1,t2,gcall,gdeepth]


 


def sim_sequencely(pMax,r):
##    ne=int(p*r)
    os.chdir("E:\\kuaipan\\python\\code")
    

    for p in range(50,pMax):
        ne=min(int(p*r),int(p*(p-1)/2))
        N=ne-p+1
        out={'n':[ ],'cout':[ ],'t_create':[ ],'t_count':[ ],'gcallout':[ ],'gdeepout':[ ]}
        
        
        G=tree_construct(p)
        global gcall
        global gdeepth

        for n in range(N):            
            for j in range(1):
                Gtemp=copy.deepcopy(G)
                t1=time.time()                
                Gtemp=perfect_non_adj_ver3(Gtemp)
                t2=time.time()
                gcall=0
                gdeepth=0
                numb=count(Gtemp,n+p)
                t3=time.time()

                out['n'].append(n+p)
                out['cout'].append(numb)
                out['t_create'].append(t2-t1)
                out['t_count'].append(t3-t2)
                out['gcallout'].append(gcall)
                out['gdeepout'].append(gdeepth)
                
            G=copy.deepcopy(Gtemp)         

        fname="chordal_p_"+str(p)+"r"+str(r)+".pkl"
        fi11=open(fname,"w")
        pickle.dump(out,fi11)
        fi11.close()
##        print out



        
def comp_two_small(p1=5,p2=14,numlab=100):
    # the number of nodes from 5 to  , edges from p(p-1)/4 to p(p-1)/2-3
    # run 100 time for each setting
    # chordal graph generated sequently
    # count is much slow when p>13 for dense chordal graph

    out=[]
    for p in range(p1,p2):        
        for numrun in range(numlab): 
            nmin=p*(p-1)/2-7
            G=random_connected_chordal_ver3(p,nmin)
            for addnum in range(1,5):
                n=nmin+addnum
                G=addEdge(G)
                t1=time.time()                
                cout1=count(G,n)
                t2=time.time()                
                t3=time.time()                
                cout2=evalsize(G)
                t4=time.time()                
                out.append([p,n,cout1,cout2,t2-t1,t4-t3])
          

    fname="com_twomethod_small"+str(p1)+"to"+str(p2)+".pkl"
    fi11=open(fname,"w")
    pickle.dump(out,fi11)
    fi11.close()
    return(out)
   
import csv
def formula_simu(p,r,numlab=100):
    # for one p and r
 
    n=p*r
    for labi in range(numlab):   
        G=random_connected_chordal_ver3(p,n)    
        t1=time.time()                
        cout1=count(G,n)
        t2=time.time()                
        cout2=evalsize(G)
        t4=time.time()           
        out=[p,r,cout2,t2-t1,t4-t2]    
        with open('comparetwo_r2to5_rep1.csv', 'ab') as fp:  
            a = csv.writer(fp, delimiter=',')
            a.writerow(out)   


def formula_simu_2(p,r,numlab=100):
    # for one p and r
    # for dominating parameter
    n=p*r
    for labi in range(numlab):   
        G=random_connected_chordal_ver3(p,n)   
        for kk in range(1,5):
            dominating=kk
            t2=time.time()                
            cout2=evalsize(G)
            t4=time.time()           
            out=[p,r,kk,cout2,t4-t2]    
            with open('comparetwo_onlyformula_r5.csv', 'ab') as fp: 
                a = csv.writer(fp, delimiter=',')
                a.writerow(out)  


    
     
def testsmallp(p,r):
    
##   change to test small p for all
    global gcall
    global gdeepth

    ne=min(int(p*r),int(p*(p-1)/2))
    N=ne-p+1
    out={'n':[ ],'cout':[ ],'t_create':[ ],'t_count':[ ]}
    
    
    Gtemp=tree_construct(p) 
    for n in range(N):        
        gcall=0
        gdeepth=0
     
        t1=time.time()                
        Gtemp=perfect_non_adj_ver3(Gtemp)
        t2=time.time() 
        numb=count(Gtemp,n+p)
        t3=time.time()

        out['n'].append(n+p)
        out['cout'].append(numb)
        out['t_create'].append(t2-t1)
        out['t_count'].append(t3-t2)
     
        print(n+p,numb,t2-t1,t3-t2)
    return out
        
                
def testlargesparse(p,n):
    
##   change to test small p for all
    global gcall
    global gdeepth

    ne=n
    N=ne-p+1
    out={'n':[ ],'cout':[ ],'t_create':[ ],'t_count':[ ]}
    gcall=0
    gdeepth=0
    
    Gtemp=tree_construct(p)
    t1=time.time()
    for n in range(N):                     
        Gtemp=perfect_non_adj_ver3(Gtemp)
    t2=time.time()
        
    print(t2-t1)  # 修复这里的 print 语法
    
    numb=count(Gtemp,n+p)
    t3=time.time()

    out['n'].append(n+p)
    out['cout'].append(numb)
    out['t_create'].append(t2-t1)
    out['t_count'].append(t3-t2)
 
    print(n+p,numb,t2-t1,t3-t2)  # 这行已经修复过了
    return out
               

  
##        print out

def mytimecomp(pMax1,pMax2,lab):
    ###p=16, sum function seems  error,gives nagtive number, but all numbers are positive
    output=[]
    os.chdir("E:\\research\\python\\test")    
    for I in range(pMax1,pMax2):
        fname="chordal_p_"+str(I)+"r"+str(5)+".pkl"
        flabel=open(fname,"r")
        out=pickle.load(flabel)
        flabel.close()
        output.append(sum(out[lab]))        
    return output



def mysimu2(p,n):
    #given p and n, show the size depth of recursion and the times of calling
    global gcall
    global gdepth
    G=tree_construct(p)
    gcall=0
    gdepth=0
    for J in range(n-p+1):        
        G=perfect_non_adj_ver3(G)
    t1=time.time()
    numb=count2(G,n,0)
    t2=time.time()
    print(t2-t1,",",gcall,gdepth)



def mysimu3(p,NT,r=100):
    #
    global gcall
    global gdepth

    ne=min(int(p*r),int(p*(p-1)/2-3))
    N=ne-p+1
    out={'n':range(N),'cout':[[] for i in range(N) ],'t_count':[[] for i in range(N) ],'depth':[[] for i in range(N)]}
    
    for J in range(NT):      
        Gtemp=tree_construct(p) 
        for n in range(N):        
            gcall=0
            gdepth=0                        
            Gtemp=perfect_non_adj_ver3(Gtemp)
            t2=time.time() 
            numb=count2(Gtemp,n+p,0)
            t3=time.time()
            out['cout'][n].append(numb)
            out['depth'][n].append(gdepth)
            out['t_count'][n].append(t3-t2)
    meanout=[]
    for key in ['cout','t_count','depth']:
         meanout.append([sum(out[key][i])/(NT*1.0) for i in range(N)])
    return meanout
    
 

