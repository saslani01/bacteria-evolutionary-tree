import blosum as bl

BLOSUM_62_MATRIX = bl.BLOSUM(62)

def protein_alignment_similarity(seq1: str, seq2: str, gap: int = -12) -> float:
    """
    This function calculates the alignment similarity score of two protein sequences using Smith Waterman algorithm.

    Args:
        seq1: protein one sequence
        seq2: protein two sequence
        gap: penalty for a gap (-12 by default)

    Returns:
        alignment_score: similarity score between the two sequences
    """

    ROW_NUMBER = len(seq1) + 1
    COL_NUMBER = len(seq2) + 1
    
    score_matrix = [[0 for _ in range(COL_NUMBER)] for _ in range(ROW_NUMBER)] # Initialize with a row and a col of zeros
    alignment_score = 0 # Keep track of max, so we do not need to find it again after filling the matrix
    
    for i in range(1, ROW_NUMBER):
        for j in range(1, COL_NUMBER):
            score = max (
                0,
                score_matrix[i - 1][j] + gap,
                score_matrix[i][j - 1] + gap,
                score_matrix[i - 1][j - 1] + BLOSUM_62_MATRIX[seq1[i - 1]][seq2[j - 1]]
            )
            
            score_matrix[i][j] = score # Update the score matrix

            alignment_score = max(score, alignment_score) # Update the best score found
       
    return alignment_score
