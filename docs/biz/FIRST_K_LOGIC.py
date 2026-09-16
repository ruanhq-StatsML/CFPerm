# manuscript first_k on trail batches
import numpy as np


def first_k_consecutive(det, k):
    """First trail index i where det[i:i+k] are all True.

    ``det`` is bool over trail batches only (after cut).
    Index 0 = first batch after cut; None = never.
    first1 = onset; first2/first3 = sustained consecutive fires.
    """
    det = np.asarray(det, dtype=bool).ravel()
    if k <= 0 or len(det) < k:
        return None
    for i in range(len(det) - k + 1):
        if det[i : i + k].all():
            return int(i)
    return None


# example
# trail_fires = [False, True, True, False]  # batches after cut
# first1 -> 1; first2 -> 1; first3 -> None
