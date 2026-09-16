def first_k_consecutive(det, k):
    det = np.asarray(det, dtype=bool)
    for i in range(len(det) - k + 1):
        if det[i:i+k].all():
            return int(i)
    return None
