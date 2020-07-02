import numpy as np
import itertools
import random	
#amplification and deletion
#with the special case of 0
def amp(seq, n):
    L = seq
    for l in range(0, np.abs(n)):
        for i in range(0, len(seq)):
            if (L[i] == 0 and n < 0):
                L[i] = 0
            else:
                if n > 0:
                    L[i] = L[i] + 1
                else:
                    L[i] = L[i] - 1
    return L

#check if amplification from zero exist 
def check_zero(start, target):
    ind = True
    for i in range(0, len(start)):
        if (start[i] == 0 and target[i] != 0):
            ind = False
    return ind

#find breakpoints for "bins"
def break_point(s1, s2):
    s1 = np.array(s1)
    s2 = np.array(s2)
    temp = s1 - s2
    bp = []
    i = 0
    for i in range(0, len(temp)-1):
        if temp[i] != temp[i+1]:
            bp.append(i+1)
        i = i + 1
    bp.insert(0, 0)
    bp.append(len(s1))
    return bp

#find all pair of start and end point from breakpoints
def find_pair(bp):
    L = len(bp)
    if (L == 2):
        return bp
    pair = []
    for j in range(1, L):
        temp1 = [bp[0], bp[j]]
        if (j != L-1):
            for k in range(0, j + 1):
                temp2 = [bp[k], bp[L-1]]
                pair.append([temp1, temp2])
                pair.append([temp2, temp1])
        else:
            for m in range(1, j):
                for n in range(1, j-m+1):
                    temp2 = [bp[m], bp[m+n]]
                    pair.append([temp1, temp2])  
                    pair.append([temp2, temp1])
    return pair




def sub_fst_max(s1, s2, bp):
    #initialization 
    dist = 0
    temp_seq = []
    k = 0
    event = []
    A1 = s1[bp[0]:bp[1]], s2[bp[0]:bp[1]]
    dif = np.abs(np.array(A1[1]) - np.array(A1[0]))
    pos = np.argmax(dif)
    dist = dif[pos]
    event = [(np.array(A1[1]) - np.array(A1[0]))[pos], bp]
    #store intermediate sequence in temp_seq
    temp_pseq = amp(A1[0], (np.array(A1[1]) - np.array(A1[0]))[pos])
    temp_seq = temp_seq + s1
    for i in range(bp[0], bp[1]):
        temp_seq[i] = temp_pseq[k]  
        k = k + 1
    #print('max', event)
    return temp_seq, dist, event 

def sub_fst_min(s1, s2, bp):
    dist = 0
    temp_seq = []
    k = 0
    #event = []
    #print(bp)
    A1 = s1[bp[0]:bp[1]], s2[bp[0]:bp[1]]
    #print('A1', A1)
    dif = np.abs(np.array(A1[1]) - np.array(A1[0]))
    #print('DIF', dif)
    pos = np.argmin(dif)
    dist = dif[pos]
    event = [(np.array(A1[1]) - np.array(A1[0]))[pos], bp]
    #print('E', event)
    #store intermediate sequence in temp_seq
    temp_pseq = amp(A1[0], (np.array(A1[1]) - np.array(A1[0]))[pos])
    temp_seq = temp_seq + s1
    for i in range(bp[0], bp[1]):
        temp_seq[i] = temp_pseq[k]  
        k = k + 1
    #print('s1', s1, 'event', event, 'temp_seq', temp_seq)
    #print('tempseq', temp_seq)
    return temp_seq, dist, event 



def sub_fst_s(s1, s2, bp):
    #print('starts')
    #print(s1, s2, bp)
    if (check_zero(s1, s2) == False):
        return s1, np.inf, []
    dist = 0
    temp_seq = []
    k = 0
    A1 = s1[bp[0]:bp[1]], s2[bp[0]:bp[1]]
    dif = np.abs(np.array(A1[1]) - np.array(A1[0]))
    #print('dif', dif)
    pos = np.argmax(dif)
    dist = dif[pos]
    event = [(np.array(A1[1]) - np.array(A1[0]))[pos], bp]
    #print('subfst', event)
    #store intermediate sequence in temp_seq
    #print(A1[0])
    temp_pseq = amp(A1[0], (np.array(A1[1]) - np.array(A1[0]))[pos])
    #print('temppseq', temp_pseq)
    temp_seq = temp_seq + s1
    #print('temp_seq', temp_seq)
    for i in range(bp[0], bp[1]):
        temp_seq[i] = temp_pseq[k]  
        k = k + 1
    return temp_seq, dist, event

#run fst 
def fst(s1, s2, bp):
    pair = find_pair(bp)
    current_dist = 0
    LCA = []
    event = []
    if (len(bp) == 2):
        final_seq, current_dist, event =  sub_fst_max(s1, s2, bp)
        LCA = [s1, s2]
        print(event)
        return LCA, current_dist, [event]
    #if 2 overlapping bin can not resolve two sequences
    if (event == []):
        for i in range(0, len(bp)-1):
            event_range = [bp[i], bp[i + 1]]
            event_op = s2[bp[i]] - s1[bp[i]]
            event.append([event_op, event_range])
            current_dist = current_dist + np.abs(event_op)
            #event = event
        #event = [event]
        LCA= [s1, s2]
    #for i in range(0, len(bp)):
        #print(bp[i]-1)
        #current_dist = current_dist + np.abs(s1[bp[i]-1]- s2[bp[i]-1])
    LCA= [s1, s2]
    for i in range(0, len(pair)):
        bp_1 = pair[i][0]
        bp_2 = pair[i][1]
        temp_seq_min, disti_1, eventi_1= sub_fst_min(s1, s2, bp_1)
        temp_seq_max, dista_1, eventa_1 = sub_fst_max(s1, s2, bp_1)
        final_seq_min, disti_2, eventi_2= sub_fst_s(temp_seq_min, s2, bp_2)
        #print(sub_fst_s(temp_seq_max, s2, bp_2))
        final_seq_max, dista_2, eventa_2= sub_fst_s(temp_seq_max, s2, bp_2)
        if (final_seq_min == s2):
            if (disti_1 + disti_2 < current_dist):
                LCA = [s1, temp_seq_min, s2]
                event = [eventi_1, eventi_2]
                current_dist = disti_1 + disti_2
            if (disti_1 + disti_2 == current_dist):
                ind = 0
                for i in range(0, len(LCA)):
                    if LCA[i] == temp_seq_min:
                        ind = 1
                if ind == 0:
                    LCA.append(temp_seq_min)
                    event.append(eventi_1)
                    event.append(eventi_2)
                    
        if (final_seq_max == s2):
            #print(s1, s2, final_seq_max, 'max', eventa_1, eventa_2)
            if (dista_1 + dista_2 < current_dist):
                LCA = [s1, temp_seq_max, s2]
                current_dist = dista_1 + dista_2
                event = [eventa_1, eventa_2]
            if (dista_1 + dista_2 == current_dist):
                ind = 0
                for i in range(0, len(LCA)):
                    if LCA[i] == temp_seq_max:
                        ind = 1
                if ind == 0:
                    LCA.append(temp_seq_max)
                    event.append(eventa_1)
                    event.append(eventa_2)
    return LCA, current_dist, event

#filter LCA contains amplification from zero
def filter(s, LCA):
    zeros = []
    L = []
    for loc in range(0, len(s)):
        if s[loc] != 0:
            zeros.append(loc)
    for i in range(0, len(LCA)):
        ind = True
        for j in range(0, len(zeros)):
            if LCA[i][zeros[j]] == 0:
                ind = False
            else:
                pass
        if ind == True:
            L.append(LCA[i])
    return L

#for general purpose only. Zero events should not exist given how seperator is defined.      
def remove_zero_event(events):
    events_updated = []
    for i in range (0, len(events)):
        if events[i][0] != 0:
                events_updated.append(events[i])
    return events_updated
            

#calculate the distance for two subsequence, the version used for infer LCA
def dist(s1, s2):
    bp = break_point(s1, s2)
    L1 = fst(s1, s2, bp)
    L2 = fst(s2, s1, bp)
    LCA = L1[0] + L2[0]
    LCA.sort()
    LCA = list(LCA for LCA,_ in itertools.groupby(LCA)) 
    LCA = filter(s1, LCA)
    LCA = filter(s2, LCA)
    return LCA, L1[1]

#calculate the distance for two subsequence, the version used for infer events
def event_dist(s1, s2):
    bp = break_point(s1, s2)
    L1 = fst(s1, s2, bp)
    L = len(s1)
    event = L1[2]
    event.sort()
    event = list(event for event,_ in itertools.groupby(event)) 
    event = remove_zero_event(event)
    return L1[1], event, L

#Compute distance and events of two sequence(has direction)
def event(s1, s2):
    LCA = np.zeros((1, len(s1)))
    distance = 0
    event = []
    temp1 = []
    temp2 = []
    if (len(s1) != len(s2)):
        print("input length does not match")
        return 0
    I = True
    for i in range (0,len(s1)):
        if s1[i]==s2[i]:
            if I == True:
                LCA[:, i] =  s1[i]
                #I = True
                pass
            else:
                #print(temp1, temp2, type(temp1))
                result = event_dist(temp1, temp2)
                distance += result[0]
                partial_event = result[1]
                r = result[2]
                for x in range(0, len(partial_event)):
                    event_op = partial_event[x][0]
                    event_temp = partial_event[x][1]
                    event_temp[0] = event_temp[0]+i-r
                    event_temp[1] = event_temp[1]+i-r
                    event.append([event_op,event_temp])
                    event_temp = []
                    event_op = []
                I = True
                distance += result[0]
                temp1 = []
                temp2 = []
        else:
            if I == True:
                I = False
                temp1.append(s1[i])
                temp2.append(s2[i])
            else:
                temp1.append(s1[i])
                temp2.append(s2[i])
    i = i + 1
    if (len(temp1) == 0):
        return LCA, distance
    result = event_dist(temp1, temp2)
    distance += result[0]
    partial_event = result[1]
    r = result[2]
    for x in range(0, len(partial_event)):
        event_op = partial_event[x][0]
        event_temp = partial_event[x][1]
        event_temp[0] = event_temp[0]+i-r
        event_temp[1] = event_temp[1]+i-r
        event.append([event_op, event_temp])
        event_temp = []
        event_op = []
    #event = list(event for event,_ in itertools.groupby(event)) 
    return distance, event        

def seq_dist(s1, s2):
    LCA = np.zeros((1, len(s1)))
    distance = 0
    temp1 = []
    temp2 = []
    if (len(s1) != len(s2)):
        print("input length does not match")
        return 0
    I = True
    for i in range (0,len(s1)):
        if s1[i]==s2[i]:
            if I == True:
                #LCA[:, i] =  s1[i]
                pass
                I = True
            else:
                result = dist(temp1, temp2)
                LCA_temp = result[0]
                if (LCA.shape[0]>16384):
                    LCA = LCA[np.random.choice(LCA.shape[0], 16384, replace=False), :]
                if (LCA.shape[0] == 16384):
                    #print(">", LCA.shape[0])
                    LCA = LCA[np.random.choice(LCA.shape[0], int(np.around(LCA.shape[0]/(len(LCA_temp)+1))), replace=False), :]
                    #print("reduce size", LCA.shape[0])
                    #sl = np.random.choice(int(LCA.shape[0]), int(np.around(LCA.shape[0]*(len(LCA_temp)-1)/len(LCA_temp), decimals =1)), replace=False)
                    #LCA= np.delete(LCA, sl, axis=0)
                current_LCA = LCA
                for j in range(0, len(LCA_temp)-1):    
                    LCA = np.append(LCA,values=current_LCA,axis=0)
                    r = len(LCA_temp[j])
                for j in range(0, len(LCA_temp)):
                    for k in range(0, r):
                        LCA[(j*current_LCA.shape[0]):((j+1)*current_LCA.shape[0]),i-r+ k]= int(LCA_temp[j][k])       
                distance += result[1]
                LCA[:, i] =  int(s1[i])
                I = True
                temp1 = []
                temp2 = []
        else:
            if I == True:
                I = False
                temp1.append(s1[i])
                temp2.append(s2[i])
            else:
                temp1.append(s1[i])
                temp2.append(s2[i])
    i = i + 1
    if (len(temp1) == 0):
        return LCA, distance
    result = dist(temp1, temp2)
    LCA_temp = result[0]

    if (LCA.shape[0]>16384):
        LCA = LCA[np.random.choice(LCA.shape[0], 16384, replace=False), :]
    if (LCA.shape[0] == 16384):
        LCA = LCA[np.random.choice(LCA.shape[0], int(np.around(LCA.shape[0]/(len(LCA_temp)+1))), replace=False), :]
        #sl = np.random.choice(int(LCA.shape[0]), int(np.around(LCA.shape[0]*(len(LCA_temp)-1)/len(LCA_temp), decimals =1)), replace=False)
        #LCA= np.delete(LCA, sl, axis=0)
    current_LCA = LCA
    for j in range(0, len(LCA_temp)-1):
        LCA = np.append(LCA,values=current_LCA,axis=0)
        r = len(LCA_temp[j])
        for j in range(0, len(LCA_temp)):
            for k in range(0, r):
                LCA[(j*current_LCA.shape[0]):((j+1)*current_LCA.shape[0]),i-r+ k]= LCA_temp[j][k]         
    distance += result[1]
    return LCA, distance

if __name__=="__main__":

    #Toy example amplification from zero 
    #The event function does not prevent amplification from zero 
    #but when infer LCAs amplification from zero is not allowed.

    X = [2,2,4,4,5,5,1,1,2,2,4,4,5,5]
    Y = [1,1,2,2,4,4,1,1,1,1,2,2,4,4]
    dists, events = event(X, Y)
    print(events)

    
