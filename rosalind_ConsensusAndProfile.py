#The profile is created by making 4 lists that correspond to a nucleotide
#One will have to count how many times a given nucleotide is called at any
#given postion.
#if a nucleotide has the most hits for a given position then it will be the
#the consensus nucleotide.
file = open('rosalind_cons.txt','r')
raw = file.readlines()
data = []
for i in raw:
    insert = i.strip('\n')
    data.append(insert)


fastanames = []
fastanamesindexes = []
dataindexes = []
for i in data:
    dataindexes.append(data.index(i))
    if '>' in i:
        fastanames.append(i)
for i in fastanames:
    index = data.index(i)
    fastanamesindexes.append(index)


dict ={}
for i in fastanamesindexes:             
    y = data[i]                         
    dict[y]={'seq':''}                                  
    dataindexes.pop(0)                  
    countforpop = 0                     
    for i in dataindexes:               
        if i not in fastanamesindexes:  
            dict[y]['seq'] = (dict[y]['seq']) + data[i]
            countforpop +=1             
            continue
        if i in fastanamesindexes:
            for i in range(0,countforpop):
                dataindexes.pop(0)      .
            break
#The three blocks of code above were used to convert the text file with a FASTA format into a dictionary
#data type. The name of the sequence is the dictionary key and the value is the DNA sequence.
#SEE my "Computing GC Content" repository for the python code that has a thorough annotation of the blocks
#of code above.


A = []                              # I made empty lists that correspond to nucleotide that is being
T = []                              #counted a specific nucleotide position; an item's index in the list
                                    #corresponds to the nucleotide postion. These lists will expand
C = []                              #until no more nucleotide postions can be read from the reference
G = []                              #sequences
for i in dict:                      #This scans through very dictionary key
    seq = dict[i]['seq']            #I want to look at one sequence at a time per dictionary key(fasta name)
    for i in range(0,len(seq)):     #The sequence is a string so I will call each character in the string
                                    #using a range because the string has repeating char and python gets
                                    #confused
        nucleotide = seq[i]
        if nucleotide == 'A':       #The list starts off as empty and I want to add a zero if a
                                    # index that corresponds to the nucleotide position I am counting in
            if len(A) < len(seq):   #for a given nucleotide,If the length of the A list is less than the length of the sequence
                                    #then that means that the nucleotide positon is not represented in the dNTP list
                                    #and I need to append a zero to all lists to synchronize the addition of
                A.append(0)         #a nucleotide postion.
                T.append(0)
                C.append(0)
                G.append(0)
                A[i] +=1
            else:
                A[i] +=1
        elif nucleotide == 'T':
            if len(T) < len(seq):   #This code could be cleaned up by defining a function that appends the 0s.
                T.append(0)
                A.append(0)
                C.append(0)
                G.append(0)
                T[i] +=1
            else:
                T[i] +=1
        elif nucleotide == 'C':
            if len(C) < len(seq):
                T.append(0)
                A.append(0)
                C.append(0)
                G.append(0)
                C[i] +=1
            else:
                C[i] +=1
        elif nucleotide == 'G':
            if len(A) < len(seq):
                T.append(0)
                A.append(0)
                C.append(0)
                G.append(0)
                G[i] +=1
            else:
                G[i] +=1



ConsensusSeq = []           #REMEMBER to put list index to the left of the equal sign when changing its value
NumList = [0,0,0,0]                             #indexes correspond from 0-3 to ATCG respectively
for i in range(0,len(A)):                         #note: make an empty list of four items--NumList. indexes in that list
    NumList[0] = A[i]           #(0-3) correspond to ATCG respectively. Those indexes will be
    NumList[1] = T[i]                  #replaced with their Nucelotide call number per column of the
    NumList[2] = C[i]      #consensus matrix which correspond to the nucleotide position
    NumList[3] = G[i]     #then find the max value in that list. trace back what index the max
    maxNum = max(NumList)                     #number corresponds to, e.g. max is at index 1 in NumList then the
    indexfornuccall = NumList.index(maxNum)                        #consensus nucleotide at that position is T
    if indexfornuccall == 0:
        ConsensusSeq.append('A')
    elif indexfornuccall == 1:
        ConsensusSeq.append('T')
    elif indexfornuccall == 2:
        ConsensusSeq.append('C')
    elif indexfornuccall == 3:
        ConsensusSeq.append('G')



print(''.join(ConsensusSeq))        #I need to convert the number list into string items before I can use 
for i in A:                         #the "join" function.
    indexref = A.index(i)               
    A[indexref] = str(i)
for i in T:
    indexref = T.index(i)
    T[indexref] = str(i)
for i in G:
    indexref = G.index(i)               
    G[indexref] = str(i)
for i in C:
    indexref = C.index(i)               
    C[indexref] = str(i)
print('A: '+ ' '.join(A)+'\n'+ 'C: ' + ' '.join(C) + '\n'+ 'G: '+' '.join(G) + '\n'+ 'T: '+' '.join(T))
