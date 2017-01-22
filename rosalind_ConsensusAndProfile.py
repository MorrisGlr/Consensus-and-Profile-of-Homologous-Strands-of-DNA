#I will be comparing multiple strands of DNA that are the same length.
#I will have to make a profile on how many times a nucleotide is called
#at a given postion in the DNA string.
#The profile is created by making 4 lists that correspond to a nucleotide
#I will have to count how many times a given nucleotide is called at any
#given postion.
#if a nucleotide has the most hits for a given position then it will be the
#the consensus nucleotide.
#it is possible that there is more than one consensus nucleotide so I will
#have to make room for that possibility.a
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
###the above block of code to identify the names of the sequences and to locate
#the names in the data list via their corresponding indexes this will be
#important for sorting between name info and actual corresponding sequencing info
#I will be manipulating the index lists while conserving the original data set list
#

dict ={}
for i in fastanamesindexes:             ###scan through all of the name indexes list
                                        #this calls the name of dict key that I want to
                                        #add seq data to for one cycle of the for loop
                                        ###
    y = data[i]                         #call the name by index in the original data set
    dict[y]={'seq':''}  #the general structure of the embedded dictionary
                                        #and it will iterate through each seq name
                                        ###
    dataindexes.pop(0)                  ###gets rid of first list item called index from
                                        #dataindexes, in this case it is the name indexes
                                        #for origanl data we do not want the first name to#
                                        #iterate in the upcoming for loop because within
                                        #that loop, name related indexes break out of the
                                        #for loop and then I do not end up inserting seq
                                        #data into the dict
                                        ###
    countforpop = 0                     #variable that is used to count the number of loops
                                        #in the for loop below it. This becomes useful in
                                        #getting rid of the data indexes that are correlated
                                        #with seqdata. Think of it as scratching off an item
                                        #in a to do list, you do not want to redo the task
    for i in dataindexes:               #iterate through the remaining data indexes that have
                                        #not been popped off yet
        if i not in fastanamesindexes:  #this checks if an index corresponds to a seq name
                                        #which is not what we want to add to the 'seq section'
                                        #for the name
                                        ###
            dict[y]['seq'] = (dict[y]['seq']) + data[i]
            countforpop +=1             #each time a seq item gets tacked onto the seq key the
                                        #the counter increases by 1 becuase I need to count how
                                        #how many items in the list I need to get rid of/ cross
                                        #off my list
            continue
        if i in fastanamesindexes:
            for i in range(0,countforpop):
                dataindexes.pop(0)      ###the for loop iterates the same number of for loop
                                        #cycles. I do this only when an index that correlates to
                                        #a name gets interating in the for loop because I want
                                        #to move onto the next name and overarching dict key
                                        #which will be the next seq name. Adding the pop function
                                        #in the 'if' condition of the for loop will mess up the
                                        #iteration and skip all around. It is hard to explain.
            break



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
            if len(T) < len(seq):
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



print(''.join(ConsensusSeq))        #I need to convert the number list into string items before I can
for i in A:
    indexref = A.index(i)               #use the "join" function
    A[indexref] = str(i)
for i in T:
    indexref = T.index(i)
    T[indexref] = str(i)
for i in G:
    indexref = G.index(i)               #use the "join" function
    G[indexref] = str(i)
for i in C:
    indexref = C.index(i)               #use the "join" function
    C[indexref] = str(i)
print('A: '+ ' '.join(A)+'\n'+ 'C: ' + ' '.join(C) + '\n'+ 'G: '+' '.join(G) + '\n'+ 'T: '+' '.join(T))
