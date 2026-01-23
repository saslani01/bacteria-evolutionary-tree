
from alignment import protein_alignment_similarity
import numpy as np

def create_similarity_matrix():
    """
    This function reads the protein_sequences.txt and creates a similarity matrix of the sequences at ./data/similarity_matrix.txt
    """
    
    with open("./data/protein_sequences.txt", 'r') as file:
        lines = file.readlines()
        sequences = [line[line.find('\t') + 1:].strip() for line in lines]

    n = len(sequences)
    matrix = [[0] * n for _ in range(n)]

    for i in range(n): 
        matrix[i][i] = protein_alignment_similarity(sequences[i], sequences[i])
        print(f"[PROCESSED BACT_{i} & BACT_{i}]")
        for j in range(i + 1, n): 
            score = protein_alignment_similarity(sequences[i], sequences[j])
            matrix[i][j] = score
            matrix[j][i] = score
            print(f"[PROCESSED BACT_{i} & BACT_{j}]")

    np.savetxt("./data/similarity_matrix.txt", np.array(matrix), fmt='%d', delimiter='\t')

def create_normalized_and_distance_matrices():
    """
    This function reads the similarity_matrix.txt and creates:
        *  normalized similarity matrix of the sequences at ./data/normalized_similarity_matrix.txt using geometric mean
        *  distance matrix at ./data/distance_matrix.txt
    """
    raw_matrix = np.loadtxt("./data/similarity_matrix.txt", delimiter = '\t')

    n = len(raw_matrix)
    normalized_matrix = [[0.0] * n for _ in range(n)]
    distance_matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        normalized_matrix[i][i] = 1.0 # diag must be 1.0
        distance_matrix[i][i] = 0.0
        print(f"[NORMALIZED AND FOUND DISTANCE BACT_{i} & BACT_{i}]")
        for j in range(i + 1, n):
            normalized_val = raw_matrix[i][j] / np.sqrt(raw_matrix[i][i] * raw_matrix[j][j]) 
            normalized_matrix[i][j] = normalized_val
            normalized_matrix[j][i] = normalized_val

            distance_matrix[i][j] = 1 - normalized_val
            distance_matrix[j][i] = 1 - normalized_val


            print(f"[NORMALIZED AND FOUND DISTANCE BACT_{i} & BACT_{j}]")
    
    np.savetxt("./data/normalized_similarity_matrix.txt", np.array(normalized_matrix), fmt='%.6f', delimiter='\t')
    np.savetxt("./data/distance_matrix.txt", np.array(distance_matrix), fmt='%.6f', delimiter='\t')

def print_furtherst_and_closest_ids():
    """
    This function calculate the furthest and closest bacterias and the distance between them.
    """
    
    distance_matrix = np.loadtxt("./data/distance_matrix.txt", delimiter = '\t')
    n = len(distance_matrix)

    min_distance = float("inf")
    max_distance = float("-inf")
    closet = None
    furthest = None
    for i in range(n):
        for j in range(i + 1, n):
            d = distance_matrix[i][j]
        
            if d < min_distance:
                min_distance = d
                closet = [i, j]
            if d > max_distance:
                max_distance = d
                furthest = [i, j]

    print(f"FURTHEST: ID's {furthest[0]} and {furthest[1]} with distance of {max_distance}\nCLOSEST: ID's {closet[0]} and {closet[1]} with distance of {min_distance}")
