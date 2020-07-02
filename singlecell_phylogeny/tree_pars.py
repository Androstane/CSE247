# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:30:34 2020

@author: eveli
"""
from ete3 import Tree
import scipy 
import ast
import numpy as np
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import average, fcluster
from distance import seq_dist

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%d%s" % (leaf_names[node.id], 0, newick)
    else:
        if len(newick) > 0:
            newick = "):%d%s" % (0, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, 0, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), 0, leaf_names)
        newick = "(%s" % (newick)
        return newick
    
def compute_tree(seq):
    L = len(seq)
    dist = np.empty((L, L))
    for i in range(0, L):
        for j in range(0, L):
            k = seq_dist(seq[i], seq[j])[1]
            if (k == np.inf):
                dist[i][j] = 999999
            else:
                dist[i][j] = k
    dist = dist[np.triu_indices(L, k = 1)]
    Z = average(dist)
    tree = hierarchy.to_tree(Z,False)
    result = getNewick(tree, "", tree.dist, seq)
    return result, tree

def score(A, B):
    if (type(A[0]) == int  and type(B[0]) == int):
        return seq_dist(A, B)
    else:
        LCA = []
        dist = np.inf
        for i in range(0, len(A)):
            for j in range(0, len(B)):
                if (seq_dist(A[i], B[j])[1] < dist):
                    dist = seq_dist(A[i], B[j])[1]
                    LCA = seq_dist(A[i], B[j])[0] 
    return LCA, dist

def find_match(p):
    p = str(p)
    s = 0
    e = 0 
    for i in range(0,len(p)):
        if (p[i] == "("):
            s = i
        if (p[i] == ")"):
            e = i
            return s, e
    return 0

def num_par(text):
    istart = []  # stack of indices of opening parentheses
    d = {}
    for i, c in enumerate(text):
        if c == '(':
             istart.append(i)
        if c == ')':
            try:
                d[istart.pop()] = i
            except IndexError:
                print('Too many closing parentheses')
    return d

def find_bracket(text):
    istart = []  # stack of indices of opening parentheses
    d = {}
    for i, c in enumerate(text):
        if c == '[':
             istart.append(i)
        if c == ']':
            try:
                d[istart.pop()] = i
            except IndexError:
                print('Too many closing parentheses')
    return d

def extract_seq(p):
    L1 = []
    L2 = []
    dist = []
    p = str(p)
    loc = find_bracket(p)
    s1 = min(loc.keys())
    e1 = loc[s1]
    e2 = max(loc.values())
    s2 = list(loc.keys())[list(loc.values()).index(e2)]
    p1 = p[s1:e1 + 1]
    p2 = p[s2:e2+1]
    temp_1 = ast.literal_eval(p1)
    if type(temp_1[0]) is int:
        L1.append(temp_1)
    else:
        L1 = temp_1
    temp_2 = ast.literal_eval(p2)
    if type(temp_2[0]) is int:
        L2.append(temp_2)
    else:
        L2 = temp_2
    for i in range(0, len(p)):
        if (p[i] == ":"):
            for j in range(i, len(p)):
                if (p[j] == "," or p[j] == ")"):
                    dist.append(float(p[i+1:j]))
                    break
    return L1, L2, dist

def pars(A,B):
    s = 0
    for i in range(0, len(A)):
        s = s + np.abs(A[i] - B[i])
    return s
#modified version of score that will return the selected LCA. 
def pars_score(A,B):
    if (type(A[0]) == int  and type(B[0]) == int):
        return seq_dist(A, B)
    else:
        LCA = []
        dist = np.inf
        ind1 = 0
        ind2 = 0
        for i in range(0, len(A)):
            for j in range(0, len(B)):
                if (seq_dist(A[i], B[j])[1] < dist):
                    dist = seq_dist(A[i], B[j])[1]
                    LCA = seq_dist(A[i], B[j])[0] 
                    ind1 = i 
                    ind2 = j
    return LCA, dist, ind1, ind2
    

def tree_pars(tree, parscore, internal_list):
    tree = str(tree)
    d = num_par(tree)
    if (len(d) == 1):
        loc = find_match(tree)
        temp = tree[loc[0]: loc[1] + 1]
        i = 0
        while (i < len(temp)):
            if (temp[i] == '\n'):
                temp = temp[0:i] + temp[i+2: len(temp)]
            else:
                i = i + 1 
        L1, L2, d = extract_seq(temp)
        #L1 is internal: 
        if (len(L1) == 2):
            root_cell = []
            for n in range(0, len(L1[0][0])):
                root_cell.append(2)
            root_cell = [root_cell]
            which = pars_score(root_cell, L1[0])[2]
            internal_list.append(L1[0][which])
            parscore = parscore + pars(L1[0][which], L1[1][0]) + pars(L1[0][which], L1[1][1]) + pars(L1[0][which], root_cell[0])
        else:
            root_cell = []
            for n in range(0, len(L1[0])):
                root_cell.append(2)
            root_cell = [root_cell]
            parscore = parscore + pars(L1[0], root_cell[0])
        #L2: 
        if (len(L2) == 2):
            which = pars_score(root_cell, L2[0])[2]
            internal_list.append(L2[0][which])
            parscore = parscore + pars(L2[0][which], L2[1][0]) + pars(L2[0][which], L2[1][1]) + pars(L2[0][which], root_cell[0])
        else:
            parscore = parscore + pars(L2[0], root_cell[0])
        return parscore
        
    else:
        loc = find_match(tree)
        temp = tree[loc[0]: loc[1] + 1]
        i = 0
        while (i < len(temp)):
            if (temp[i] == '\n'):
                temp = temp[0:i] + temp[i+2: len(temp)]
            else:
                i = i + 1 
        L1, L2, d = extract_seq(temp)  
        #both are internal nodes 
        if (len(L1) == 2 and len(L2) == 2):
            LCA, dist, t1, t2 = pars_score(L1[0], L2[0])
            parscore = parscore + pars(L1[0][t1], L1[1][0]) + pars(L1[0][t1], L1[1][1]) + pars(L2[0][t2], L2[1][0])+ (pars(L2[0][t2], L2[1][1]))
            LCA = [l.astype(int).tolist() for l in LCA]
            LCA = [LCA, [L1[0][t1], L2[0][t2]]]
            internal_list.append(L1[0][t1])
            internal_list.append(L2[0][t2])
        #L1 is internal
        elif (len(L1) == 2):
            LCA, dist, t1, t2 = pars_score(L1[0], L2)
            parscore = parscore + pars(L1[0][t1], L1[1][0]) + pars(L1[0][t1], L1[1][1])
            LCA = [l.astype(int).tolist() for l in LCA]
            LCA = [LCA, [L1[0][t1], L2[t2]]]
            internal_list.append(L1[0][t1])
        #L2 is internal
        elif (len(L2) == 2):
            LCA, dist, t1, t2 = pars_score(L1, L2[0])
            parscore = parscore + pars(L2[0][t2], L2[1][0]) + pars(L2[0][t2], L2[1][1])
            LCA = [l.astype(int).tolist() for l in LCA]
            LCA = [LCA, [L1[t1], L2[0][t2]]]
            internal_list.append(L2[0][t2])
        else:
            LCA, dist, t1, t2 = pars_score(L1, L2)
            LCA = [l.astype(int).tolist() for l in LCA]
            LCA = [LCA, [L1[t1], L2[t2]]]
        tree = str(tree[0:loc[0]] + str(LCA) + ":" + str(dist + d[0] + d[1]) + tree[loc[1] + 3 : len(tree)])
        return tree_pars(tree, parscore, internal_list)
                    
                    
def tree_score(tree):
    tree = str(tree)
    d = num_par(tree)
    if (len(d) == 1):
        loc = find_match(tree)
        temp = tree[loc[0]: loc[1] + 1]
        i = 0
        while (i < len(temp)):
            if (temp[i] == '\n'):
                temp = temp[0:i] + temp[i+2: len(temp)]
            else:
                i = i + 1 
        L1, L2, d = extract_seq(temp)
        root_cell = []
        for n in range(0, len(L1[0])):
            root_cell.append(2)
        root_cell = [root_cell]
        return d[0] + d[1] + score(L1, root_cell)[1] + score(L2, root_cell)[1]
        
    else:
        loc = find_match(tree)
        temp = tree[loc[0]: loc[1] + 1]
        i = 0
        while (i < len(temp)):
            if (temp[i] == '\n'):
                temp = temp[0:i] + temp[i+2: len(temp)]
            else:
                i = i + 1 
        L1, L2, d = extract_seq(temp)  
        LCA, dist = score(L1, L2)
        LCA = [l.astype(int).tolist() for l in LCA]
        tree = str(tree[0:loc[0]] + str(LCA) + ":" + str(dist + d[0] + d[1]) + tree[loc[1] + 3 : len(tree)]) 
        return tree_score(tree)

def infer_tree(tree):
    ilist = []
    s = 0
    score = tree_pars(tree, s, ilist)
    return score, ilist

#add [] to name of the node to newick tree string 
def add_bracket(newick_tree):
    L = len(newick_tree)
    i = 0
    ind = False
    s = newick_tree
    while i < L:
        #indicator for adding right bracket 
        #insert left bracket 
        if s[i] == 'l' or s[i] == 'D':
            s = s[0:i] + "[" + s[i:L]
            L = L+1
            i = i+2
            ind = True
        elif s[i] == ':' and ind == True:
            s = s[0:i] + "]"+s[i:L]
            L = L + 1
            i = i + 2
            ind = False
            
        else:
            i = i + 1
    return s
        
#return the order of internal node inferred from infer_tree
def name_order(tree, order):
    d = num_par(tree)
    if (len(d) == 1):
        loc = find_match(tree)
        temp = tree[loc[0]: loc[1] + 1]
        brac = find_bracket(temp)
        c_node = ""
        for key in brac:
            c_node = c_node + temp[int(key)+1: int(brac[key])]
        order.append(c_node)
        return order
    else:
        loc = find_match(tree)
        temp = tree[loc[0]: loc[1] + 1]
        brac = find_bracket(temp)
        c_node = ""
        for key in brac:
            c_node = c_node + temp[int(key)+1: int(brac[key])]
        order.append(c_node)
        tree = tree[0:loc[0]] + '[' +  c_node + ']' + tree[loc[1]+1 : len(tree)]
        return name_order(tree, order)
    
def node_order(tree):
    tree = add_bracket(tree)
    node = []
    name_order(tree, node)
    return node
if __name__=="__main__":

    #Toy example 
    X = [[2,1,1,1,2],[4,5,4,4,3],[3,0,2,5,3],[1,0,0,4,3],[2,3,3,4,2],[0,1,1,2,1]]
    Y = "([2, 1, 1, 1, 2]:0,(([0, 1, 1, 2, 1]:0,[2, 3, 3, 4, 2]:0):0,(([1, 0, 0, 4, 3]:0,[3, 0, 2, 5, 3]:0):0,[4, 5, 4, 4, 3]:0):0):0)"
    Z =  "((([4, 5, 4, 4, 3]:0,([3, 0, 2, 5, 3]:0,[1, 0, 0, 4, 3]:0):0):0,([2, 3, 3, 4, 2]:0,[0, 1, 1, 2, 1]:0):0):0,[2, 1, 1, 1, 2]:0)"
    S = "(((leaf1:0,(leaf2:0,leaf4:0):0):0,(leaf5:0,leaf6:0):0):0,D:0)"
    # X = [[1, 1, 1, 6, 6], [1, 0, 1, 3, 3], [1, 0, 0, 1, 2], [1, 3, 5, 5, 3]]
    #string, tr = compute_tree(X)
    #print(string)
    #print(tr)
    #print(tree_score(string))
    #print("tree1")
    print(infer_tree(Y))
    print(node_order(S))
    #print("tree2")
    #print(tree_score(Y))
