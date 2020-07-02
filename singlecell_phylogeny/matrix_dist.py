import numpy as np
def diff(a, b):
    if a > b: 
        return 1
    if a < b: 
        return -1
    else: 
        return 0

def pairwise_dist(a):
    l = a.shape[0]
    d = np.empty([l, l])
    for i in range(0, l):
        for j in range(0, l):
            print(a[i], a[j])
            d[i][j] = diff(a[i],a[j])
    return d
            
        
#Compute the percentage of the relative relationship (dist(A,B)>dist(A,C)) between leaves captured corrrectly     
def rel(a, b):
    diff = 0
    a = a.flatten()
    b = b.flatten()
    print(a, b)
    a = pairwise_dist(a)
    b = pairwise_dist(b)
    L = a.shape[0]
    for i in range(0, L):
        for j in range(0, L):
            if a[i][j] == b[i][j]:
                pass
            else:
                diff = diff + 1
    return 1 - diff/(L*L)

def fb(a,b):
    L = a.shape[0]
    s = 0
    for i in range(0, L):
        for j in range(0, L):
            s = s + np.square(a[i][j] - b[i][j])
    return s




m1 = np.array([[2,3], [2,2]])
m2 = np.array([[2,4], [2,2]])
print(rel(m1,m2))







