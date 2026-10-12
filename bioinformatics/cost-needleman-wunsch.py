from collections.abc import Callable


def cost_needleman_wunsch(seq_1: str, seq_2: str, score_fn: Callable[[str, str], int]) -> float:
    """Calculate optimal alignement score using cost-only Needleman-Wunsch algorithm
    Time complexity: O(nm)
    Space complexity: O(m)

    :param seq_1: The first sequence
    :type seq_1: str
    :param seq_2: The second sequence
    :type seq_2: str
    :param score_fn: Score function used to define the score for a match, mismatch and indel for 2 characters being compared
    :type score_fn: Callable[[str, str], int]
    :returns: The alignment score
    :rtype: float
    """
    n, m = len(seq_1), len(seq_2)
    values = [[ float("-inf") for _ in range(m + 1)] for _ in range(2)] #O(m) space, keeping 2 rows

    #init table
    values[0][0] = 0
    for c in range(1, m + 1):
        values[0][c] = values[0][c - 1] + score_fn("_", seq_2[c - 1])

    # iteratively generate alignment scores
    for r in range(1, n + 1):
        values[1][0] = values[0][0] + score_fn(seq_1[r - 1], "_")
        for c in range(1, m + 1):
            left: float = max(values[1][c], values[1][c - 1] + score_fn("_", seq_2[c - 1])) # left case
            up: float = max(values[1][c], values[0][c] + score_fn(seq_1[r - 1], "_")) # up case
            diag: float = max(values[1][c], values[0][c - 1] + score_fn(seq_1[r - 1], seq_2[c - 1])) # diagonal case
            values[1][c] = max(left, up, diag)
        values[0] = list(values[1])
        values[1] = [float("-inf") for _ in range(m + 1)]
    return values[0][m]

def score_fn(char_1: str, char_2: str) -> int:
    if char_1 == "_" or char_2 == "_":
        # indel case
        return -1
    if char_1 == char_2:
        # match case
        return 2
    else:
        # mismatch case
        return -1

if __name__ == "__main__":
    seq_1: str = "TT"
    seq_2: str = "AA"
    score: float = cost_needleman_wunsch(seq_1, seq_2, score_fn)
    print(f"Alignment score: {score}")

