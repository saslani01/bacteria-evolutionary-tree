from protein_seq_loader import create_similarity_matrix, create_normalized_and_distance_matrices
from visual import visualize

def main():
    create_similarity_matrix()   
    create_normalized_and_distance_matrices()
    visualize()

if __name__ == "__main__":
    main()
    