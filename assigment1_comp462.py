# Sequence alignment slippage-aware program
# by Theodoro Gasperin Terra Camargo 260842764

import sys
sys.setrecursionlimit(1500) # Extend recursion limit of system

#------------------------------ Classes ----------------------------------------

class Graph:

    def __init__(self, S, T):
        self.root = None
        self.node_n = 0
        self.seq1 = S
        self.seq2 = T
        self.list_of_nodes = [] 

    def add_node(self, n_parent, i, j, score):

        if self.root == None:
            self.root = Node(n_parent, 0, 0, 0)
            self.list_of_nodes.insert(self.node_n, self.root)
            self.node_n += 1

        else:
            new_node = Node(n_parent, i, j, score)
            self.list_of_nodes.insert(self.node_n, new_node)
            self.node_n += 1
    

class Node:

    def __init__(self, n_parent, i, j, score):
        self.n_parent = n_parent    # Parent index in list of Nodes in the Graph object
        self.i = i                  # Index i at Xij
        self.j = j                  # Index j at Xij
        self.score = score          # Best allignment score
    
    def print_node(self):
        print("Node: " + str(self.i) + " , " + str(self.j) + " , " + str(self.score))


# ------------------------------ Methods -----------------------------------------

def best_alignment_matrix(S, T, match, missmatch, cs, cn):
    '''
    S => Sequence of nucleotides; T => Sequence of Nucleotides
    Initializes a best alignment matrix with all values of 0, given S and T. 
    Also creates Graph of nodes used compute the optimal alignment of S vs T
    '''
    
    # Creating graph object: 
    myGraph = Graph(S, T) 
    # Xij Matrix of alignment scores:
    grid = []   

    i = 0
    while i < len(S) + 1:

        grid_line = []
        
        j = 0 
        while j < len(T) + 1:  
                     
            # First line and first collum cases
            if i == 0 and j == 0:
                grid_line.append(0)
                myGraph.add_node(None, 0, 0, 0) 
            
            elif i == 0 and j > 0:
                if T[j-1-1] == T[j -1]:# Slipage Case 1: T sequence has repeated bases Add "-" to S alignment
                    best_score = grid_line[j-1] + cs
                else:
                    best_score = grid_line[j-1] + cn
                grid_line.append(best_score)
                myGraph.add_node(i*(len(T)+1)+ (j-1), i, j, best_score)

            elif i > 0 and j == 0:
                if S[i-1-1] == S[i -1]:# Slipage Case 2: S sequence has repeated bases Add "-" to T alignment
                    best_score = grid[i-1][j] + cs
                else:
                    best_score = grid[i-1][j] + cn
                
                grid_line.append(best_score)
                myGraph.add_node(((i-1)*(len(T)+1)+ j), i, j, best_score)

            # All other cases
            elif i > 0 and j > 0 :
                #case1 -> diagonal 
                #case2 -> up                                   
                #case3 -> left

                case1 = grid[i-1][j-1] + calculate_score(S, T, i -1 , j -1, match, missmatch) # i - 1 and j -1 bcs I added a zero in front of the strings
                                                
                if S[i-1-1] == S[i -1]:# Slipage Case 2: S sequence has repeated bases Add "-" to T alignment
                    case2 = grid[i-1][j] + cs
                else:
                    case2 = grid[i-1][j] + cn
                if T[j-1-1] == T[j -1]:          # Slipage Case 1: T sequence has repeated bases Add "-" to S alignment
                    case3 = grid_line[j-1] + cs  # index i out of bounds for grid list bcs line not added yet to grid therefore grid_line used here
                else:
                    case3 = grid_line[j-1] + cn

                best_score = max(case1, max(case2, case3))  # max out of the 3 cases to add to matrix cell
                grid_line.append(best_score)
                
                # Adding Nodes to Graph:
                if best_score == case1:
                    parent_number = (i-1)*(len(T)+1) + j-1 # Calculates parent number node 
                    myGraph.add_node(parent_number, i, j, best_score) # Diagonal
                    
                elif best_score == case2:
                    parent_number = (((i-1)*(len(T)+1)) + j)
                    myGraph.add_node(parent_number , i, j, best_score) # Up
                    
                elif best_score == case3:
                    parent_number =    i*(len(T) + 1) + j - 1 
                    myGraph.add_node(parent_number, i, j, best_score) # Left             
                    
            j +=1
        grid.append(grid_line)  # Appending line to matrix
        i +=1
    
    
    aligment = print_alignment(myGraph.list_of_nodes, myGraph.list_of_nodes[-1], S, T, "", "")
    #print_matrix(grid)

    print("\nOptimal Alignment Score: " + str(grid[len(grid)-1][len(grid[0])-1]) )

    return "\nGlobal Alignment:\n" + aligment


def print_alignment(node_list, node, S, T, aligment1, aligment2):
    '''
    Given list of nodes, Strings S and T, 
    This function computes the global alignment of seq S and T
    '''

    if node.n_parent == None:
        result = "]" + aligment2 + "[\n]" + aligment1 + "["
        return result[::-1]
    else:
        next_node = node_list[node.n_parent]

        if next_node.i ==  node.i - 1 and next_node.j == node.j - 1: # Diagonal -> add nucleotide to both alignments 
            aligment1 = aligment1 + S[node.i-1]
            aligment2 = aligment2 + T[node.j-1]

        if next_node.i ==  node.i - 1 and next_node.j == node.j:     # Up   -> add nucleotide to alignment1 and gap "-" to alig2
            aligment1 = aligment1 + S[node.i-1]                         
            aligment2 = aligment2 + "-"

        if next_node.i ==  node.i and next_node.j == node.j - 1:     # Left -> add nucleotide to alignment2 and gap "-" to alig1
            aligment1 = aligment1 + "-"
            aligment2 = aligment2 + T[node.j-1]
            
        return print_alignment(node_list, next_node, S, T, aligment1, aligment2)


def calculate_score(S, T, i, j, match, missmatch):
    '''Calculates pair Alignment Score of chars in S and T, given index of both strings'''
    
    #ACTG VS ACTG (M)
    #substitution_cost_matrix = [ 
    #            [ 1, -1, -1, -1],
    #            [-1,  1, -1, -1], 
    #            [-1, -1,  1, -1],
    #            [-1, -1, -1,  1]] 

    #if S[i] == "A":
    #    l = 0
    #elif S[i] == "C":
    #    l = 1
    #elif S[i] == "T":
    #    l = 2
    #elif S[i] == "G":
    #    l = 3

    #if T[j] == "A":
    #    p = 0
    #elif T[j] == "C":
    #    p = 1
    #elif T[j] == "T":
    #    p = 2
    #elif T[j] == "G":
    #    p = 3

    #score = substitution_cost_matrix[l][p]

    if S[i] == T[j]:
        score = match
    else:
        score = missmatch

    return score

def print_matrix(matrix):
    '''Prints Matrix of Alignment scores'''
    
    print("Matrix:")
    for line in matrix:
        print("[ ", end="")
        for cell in line:
            print( str(cell) + ", ", end="")
        print("]\n", end="")

def compute_global_alignment(fasta_file, match_score, missmatches, slip_gap_penalty, non_slip_gap_penalty):

    fobj = open(fasta_file, "r", encoding="utf-8")
    sequences = []
    n_seq = 0
    
    # Extracting Nucleic Acid Sequences
    for line in fobj:

        if line == "\n":
            continue
        if ">" in line:
            sequences.insert(n_seq,"")
            continue  
        else:
            seq = sequences[n_seq]
            for i in line:
                if i == "\n":
                    continue
                seq = seq + i
            
            sequences[n_seq] = seq  #Updating entry
        
        n_seq += 1
    
    print(best_alignment_matrix(sequences[0], sequences[1], match_score, missmatches, slip_gap_penalty, non_slip_gap_penalty))



#--------------------- Main --------------------------
if __name__ == '__main__':
    #seq2 = "CAPE"
    #seq1 = "APPLE"
    #print(best_alignment_matrix(seq1, seq2, 1, -1, -1, -2))
    compute_global_alignment("big_seq.fa", 1, -1, -1, -2)




