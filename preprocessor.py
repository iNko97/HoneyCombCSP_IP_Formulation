import numpy as np

# Factory settings

L_min = 700  # Minimum panel length
L_max = 3100  # Maximum panel length
n_s_max = 2  # Maximum number of stock sizes
n_w_max = 2  # Maximum number of widths
W = [1200, 1400, 1550, 1600]  # Set of available stock widths

# Order Data: Define item types with their Width, Length, and Demand
order_data = [
    [230, 250, 600],
    [230, 2140, 600],
    [290, 2140, 600]
]


def powerset_index(index, original_set):

    binary_repr = format(index, f'0{len(original_set)}b')

    subset = [original_set[i] for i in range(len(original_set)) if binary_repr[i] == '1']

    return subset


def find_indexes_with_element(element, original_set):
    indexes = []

    if element not in original_set:
        return []

    element_idx = original_set.index(element)
    n = len(original_set)

    for combination in range(2 ** (n - 1)):

        binary_repr = format(combination, f'0{n - 1}b')

        new_binary_repr = binary_repr[:element_idx] + '1' + binary_repr[element_idx:]


        indexes.append(int(new_binary_repr, 2))

    return indexes


# 3. Define potential stocks (J)
J = np.zeros((len(W), n_s_max))

# 4. Define subsets of I compatible with stock size j (C_j)
C_j = {j: np.empty for j in J}

for j in J:
    width_j = j

    for subset in I_c:
        subset_indices = np.isin(item_names, subset)
        if np.all(item_widths[subset_indices] <= width_j):
            C_j[j].append(subset)

# 5. Define dictionary parameter a_ic indicating if item i is in subset c (a_ic)
a_ic = {(i, c): 1 if i in c else 0 for c in I_c for i in I}

# 6. Define minimum length of stock lmin_cjk for satisfying demand (lmin_cjk)
lmin_cjk = {}
for j in J:
    width_j, length_j = j
    for c in C_j[j]:
        for k in range(1, 10):  # Arbitrary upper limit for k, can be tuned
            # Get indices for the current subset `c`
            subset_indices = np.isin(item_names, c)
            # Compute maximum length required for this subset `c` and `k` panels
            demands = item_demands[subset_indices]
            lengths = item_lengths[subset_indices]
            max_len = np.max(lengths * np.ceil(demands / k))
            lmin_cjk[(c, j, k)] = max(700, min(max_len, 3100))

# 7. Define bounds for number of panels (kmin_cj, kmax_cj)
kmin_cj = {}
kmax_cj = {}
for j in J:
    for c in C_j[j]:
        # Get indices for the current subset `c`
        subset_indices = np.isin(item_names, c)
        lengths = item_lengths[subset_indices]
        demands = item_demands[subset_indices]

        # Calculate cmax and cmin values
        cmax_i = np.floor(3100 / lengths)
        cmin_i = np.floor(np.maximum(700, lengths) / lengths)

        # Calculate kmin and kmax for each subset `c` and stock size `j`
        kmin_ijn_i = np.ceil(demands / cmax_i)
        kmax_ijn_i = np.ceil(demands / cmin_i)

        kmin_cj[(c, j)] = int(np.max(kmin_ijn_i))
        kmax_cj[(c, j)] = int(np.max(kmax_ijn_i))

# Print results for validation
print("Set of item types (I):", I)
print("Subsets of item types (I_c):", I_c)
print("Potential stocks (J):", J)
print("Stocks with specific widths (J_w):", J_w)
print("Subsets of I compatible with stock size j (C_j):", C_j)
print("Parameter a_ic:", a_ic)
print("Minimum length of stock lmin_cjk:", lmin_cjk)
print("Bounds for number of panels (kmin_cj, kmax_cj):", kmin_cj, kmax_cj)


