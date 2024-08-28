import numpy as np
from itertools import product

# Factory settings

L_min = 700  # Minimum panel length
L_max = 3100  # Maximum panel length
n_s_max = 2  # Maximum number of stock sizes
n_w_max = 2  # Maximum number of widths
W = [1200, 1400, 1550, 1600]  # Set of available stock widths

# I: item types with their Width, Length, and Demand
I = [
    [230, 250, 600],
    [230, 2140, 600],
    [290, 2140, 600]
]

#Returns the I subset corresponding to a specified index.
def powerset_at_index(index):

    binary_repr = format(index, f'0{len(I)}b')

    subset = [I[i] for i in range(len(I)) if binary_repr[i] == '1']

    return subset

#Returns all the indexes that contain a certain I element
def find_indexes_with_element(element):
    indexes = []

    if element not in I:
        return []

    element_idx = I.index(element)
    n = len(I)

    for combination in range(2 ** (n - 1)):

        binary_repr = format(combination, f'0{n - 1}b')

        new_binary_repr = binary_repr[:element_idx] + '1' + binary_repr[element_idx:]

        indexes.append(int(new_binary_repr, 2))

    return indexes

#Returns all possible non-dominated N_c given an i in I_c and a w in W
def n_c_generator(_index, _w_max):

    items = powerset_at_index(_index)
    widths = [item[0] for item in items]

    max_rows = [_w_max // _width for _width in widths]

    possible_combinations = product(*(range(1, mr + 1) for mr in max_rows))

    valid_combinations = []
    for combination in possible_combinations:
        _total_width = sum(widths[i] * combination[i] for i in range(len(widths)))

        if _total_width <= _w_max:
            # Check if the new combination is dominated
            if not is_dominated(combination, valid_combinations):
                valid_combinations = remove_dominated(valid_combinations, combination)
                valid_combinations.append(combination)

    return valid_combinations

def is_dominated(new_combo, valid_combinations):
    for combination in valid_combinations:
        # Check if combo dominates new_combo
        if all(combination[i] >= new_combo[i] for i in range(len(new_combo))):
            return True
    return False

def remove_dominated(valid_combinations, new_combo):
    return [combination for combination in valid_combinations if not all(new_combo[i] >= combination[i] for i in range(len(combination)))]


# 3. Define potential stocks (J) of dimensions W * n_s_max
# each row indexes are ordered by available widths.
J = np.zeros((len(W), n_s_max))

# 4. Define subsets of I compatible with stock size j (C_j)
# Although referring to the J matrix, it only contains W values.
# J[idx, _] --> C_j[idx]
C_j = [[] for _ in range(len(W))]

for idx, width in enumerate(W):
    for index in range(1, 2**len(I)+1):
        subset = powerset_at_index(index)
        total_width = sum(item[0] for item in subset)

        if total_width <= width:
            C_j[idx].append([index])



#TODO THEN compute lmin_cjk foreach N_c, THEN get minimum



# Print results for validation
print("Set of item types (I):", I)
print("Subsets of item types (I_c):", I_c)
print("Potential stocks (J):", J)
print("Stocks with specific widths (J_w):", J_w)
print("Subsets of I compatible with stock size j (C_j):", C_j)
print("Parameter a_ic:", a_ic)
print("Minimum length of stock lmin_cjk:", lmin_cjk)
print("Bounds for number of panels (kmin_cj, kmax_cj):", kmin_cj, kmax_cj)

