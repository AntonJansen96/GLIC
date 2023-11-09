import copy

def makeSuperDict(keyLists: list):
    """Initializes a (recursively) nested dictionary of dictionaries based on a list
    of lists of keys. Helpful for organizing large amounts of data.

    Example:
        keyLists = [['A', 'B'], [1, 2], 0.33] -> {'A': {1: 0.33, 2: 0.33}, 'B': {1: 0.33, 2: 0.33}}.

    Args:
        keyLists (list): the list of lists of nested keys. The last element corresponds
        to the value of the innermost key-value pair.

    Returns:
        dict: desired nested dictionary or 'superDict' structure.
    """

    # Assertions
    assert isinstance(keyLists[0], list)
    assert len(keyLists) > 1

    # Calculate size of the superDict and provide user update.
    size = 1
    for idx in range(0, len(keyLists) - 1):
        size *= len(keyLists[idx])
    print('Created superDict of size {} (recursion depth {}).'.format(size, len(keyLists) - 1))

    # EXPLANATION THROUGH THE FOLLOWING EXAMPLE:

    # A = []
    # B = {}
    # for key in metrics:
    #     B[key] = copy.deepcopy(A)
    # C = {}
    # for key in chains:
    #     C[key] = copy.deepcopy(B)
    # D = {}
    # for key in reps:
    #     D[key] = copy.deepcopy(C)
    # E = {}
    # for key in sims:
    #     E[key] = copy.deepcopy(D)

    # List  = [sims, reps, chains, metrics, []]
    # array = [{}, {}, {}, {}, List[4]]

    # for key in List[3]:
    #     array[3][key] = copy.deepcopy(array[4])

    # for key in List[2]:
    #     array[2][key] = copy.deepcopy(array[3])

    # for key in List[1]:
    #     array[1][key] = copy.deepcopy(array[2])

    # for key in List[0]:
    #     array[0][key] = copy.deepcopy(array[1])

    array = [copy.deepcopy({}) for _ in range(0, len(keyLists) - 1)]
    array += [copy.deepcopy(keyLists[-1])]

    for idx in range(1, len(keyLists))[::-1]:
        for key in keyLists[idx - 1]:
            array[idx - 1][key] = copy.deepcopy(array[idx])

    return array[0]
