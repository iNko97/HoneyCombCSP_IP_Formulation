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


# Returns the I subset corresponding to a specified index.
def powerset_at_index(index):
    binary_repr = format(index, f'0{len(I)}b')

    subset = [I[i] for i in range(len(I)) if binary_repr[i] == '1']

    return subset


def n_c_generator(items, _w_max):
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
    return [combination for combination in valid_combinations if
            not all(new_combo[i] >= combination[i] for i in range(len(combination)))]


def function_name(index):
    for _width in W:
        subset = powerset_at_index(index)
        _demands = [item[2] for item in subset]


        for n_c in n_c_generator(subset, _width):
            # ceiling(a/b) = -floor(a/-b)
            # c_{ijn_i}
            columns_to_satisfy_demand = [-(_demands[i] // -n_c[i]) for i in range(len(_demands))]
            # cmax_{i}
            maximum_columns = [L_max // item[1] for item in subset]
            # kmin_{ijn_i}
            minimum_stock_sizes = [-(columns_to_satisfy_demand[i] // -maximum_columns[i])
                                   for i in range(len(columns_to_satisfy_demand))]
            # cmin_{i}
            minimum_columns = [max(L_min, item[1]) // item[1] for item in subset]
            # kmax_{ijn_i}
            maximum_stock_sizes = [-(columns_to_satisfy_demand[i] // -minimum_columns[i])
                                   for i in range(len(columns_to_satisfy_demand))]

            # kmin_cj
            minimum_stock_size = max(minimum_stock_sizes)
            # kmax_cj
            maximum_stock_size = max(maximum_stock_sizes)

            if not best_minimum_stock_size_length.size:
                best_minimum_stock_size_length = np.full(maximum_stock_size + 1, np.inf)
                print("Creating array!")

            for k in range(minimum_stock_size, maximum_stock_size + 1):
                # chat_{ijn_ik}
                columns_per_panel = [-(column_to_satisfy_item_demand // -k)
                                     for column_to_satisfy_item_demand in columns_to_satisfy_demand]
                item_minimum_stock_size_length = [subset[i][1] * columns_per_panel[i] for i in range(len(subset))]
                maximum_item_minimum_stock_size_length = max(item_minimum_stock_size_length)
                print(maximum_item_minimum_stock_size_length)
                print("current best", best_minimum_stock_size_length[k])
                if maximum_item_minimum_stock_size_length > L_max:
                    print(
                        f"Infeasible to meet the demand of I_{index} with {k} panels of width {_width} with row distribution of {n_c}")
                else:
                    if maximum_item_minimum_stock_size_length < best_minimum_stock_size_length[k] and maximum_item_minimum_stock_size_length > L_min:
                        np.put(best_minimum_stock_size_length, k, maximum_item_minimum_stock_size_length)
                        print("new best", type(maximum_item_minimum_stock_size_length),
                              maximum_item_minimum_stock_size_length)


function_name(8)
