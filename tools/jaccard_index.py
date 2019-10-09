def jaccard_index(v1, v2):
    intersection = []
    union = []
    for i in range(len(v1)):
        a, b = v1[i], v2[i]
        intersection.append(min(a, b))
        union.append(max(a, b))
    intersection_sum = sum(intersection)
    union_sum = sum(union)
    j = intersection_sum / union_sum
    return j