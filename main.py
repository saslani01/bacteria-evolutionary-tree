from protein_seq_loader import create_similarity_matrix, create_normalized_and_distance_matrices, print_furtherst_and_closest_ids
from visual import visualize

def main():
    create_similarity_matrix()   
    create_normalized_and_distance_matrices()
    visualize()
    print_furtherst_and_closest_ids()

if __name__ == "__main__":
    main()
    