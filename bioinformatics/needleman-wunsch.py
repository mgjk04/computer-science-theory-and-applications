from _collections_abc import Callable


def needleman_wunsch(seq_1: str, seq_2: str, score_fn: Callable[[str, str], int]) -> tuple[float, str, str]:
    """Returns optimal global alignment using the needleman-wunsch algorithm

    Time complexity: O(nm)

    Space complexity: O(nm)

    :param seq_1: The first sequence
    :type seq_1: str
    :param seq_2: The second sequence
    :type seq_2: str
    :param score_fn: Score function used to define the score for a match, mismatch and indel for 2 characters being compared
    :type score_fn: Callable[[str, str], int]
    :returns: A tuple of the alignment score and a pair of possible aligments
    :rtype: tuple[float, str, str]

    """
    n, m= len(seq_1), len(seq_2)
    values: list[list[float]]= [[float("-inf") for _ in range(m + 1)] for _ in range(n + 1)] # n+1 x m+1
    path: list[list[tuple[int, int]]] = [[(0, 0) for _ in range(m + 1)] for _ in range(n + 1)] # (0, 0) no dir, else same convention as dir
    dirs: list[tuple[int, int]] = [(-1, 0), (0, -1), (-1, -1)] # dir up, left & diagonal

    # init table
    values[0][0] = 0
    for r in range(1, n + 1):
        values[r][0] = values[r - 1][0] + score_fn(seq_1[r - 1], "_")
        path[r][0] = dirs[0]
    for c in range(1, m + 1):
        values[0][c] = values[0][c - 1] + score_fn("_", seq_2[c - 1])
        path[0][c] = dirs[1]

    # populate table
    for r in range(1, n + 1):
        for c in range(1, m + 1):
            if values[r - 1][c] + score_fn(seq_1[r - 1], "_") > values[r][c]:
                #up case
                values[r][c] = values[r - 1][c] + score_fn(seq_1[r - 1], "_")
                path[r][c] = dirs[0]
            if values[r][c - 1] + score_fn("_", seq_2[c - 1]) > values[r][c]:
                #left case
                values[r][c] = values[r][c - 1] + score_fn("_", seq_2[c - 1])
                path[r][c] = dirs[1]
            if values[r - 1][c - 1] + score_fn(seq_1[r - 1], seq_2[c - 1]) > values[r][c]:
                #diagonal case
                values[r][c] = values[r - 1][c - 1] + score_fn(seq_1[r - 1], seq_2[c - 1])
                path[r][c] = dirs[2]

    # backtrack
    align_1, align_2 = "", ""
    r, c = n, m
    while r > 0 or c > 0:
        d: tuple[int, int] = path[r][c]
        if d == dirs[0]:
            # up case
            align_1 += seq_1[r - 1]
            align_2 += "_"
        if d == dirs[1]:
            align_2 += seq_2[c - 1]
            align_1 += "_"
            # left case
        if d == dirs[2]:
            # diag case
            align_1 += seq_1[r - 1]
            align_2 += seq_2[c - 1]
        r += d[0]
        c += d[1]
    return values[n][m], align_1[::-1], align_2[::-1]

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
    seq_1: str = "ACAATCC"
    seq_2: str = "AGCATGC"
    score, align_1, align_2 = needleman_wunsch(seq_1, seq_2, score_fn)
    print(f"Alignment score: {score}\nAlignment 1: {align_1}\nAlignment 2: {align_2}")
