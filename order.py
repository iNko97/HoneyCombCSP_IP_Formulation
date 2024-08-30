import numpy as np
from itertools import product
import pandas as pd

# Factory settings
#TODO: IMPORT FROM CSV

def is_dominated(new_combo, valid_combinations):
    for combination in valid_combinations:
        # Check if combo dominates new_combo
        if all(combination[i] >= new_combo[i] for i in range(len(new_combo))):
            return True
    return False


def remove_dominated(valid_combinations, new_combo):
    return [combination for combination in valid_combinations
            if not all(new_combo[i] >= combination[i] for i in range(len(combination)))]

# Checks whether item i is in index c
def a_ic_generator(i, c):
    return 1 if bool(c & (1 << i)) else 0


class Order:
    

    # Returns the I subset corresponding to a specified index.
    def powerset_at_index(self, index):

        binary_repr = format(index, f'0{len(self.Items)}b')

        subset = [self.Items[i] for i in range(len(self.Items)) if binary_repr[i] == '1']

        return subset

    # Returns all the indexes that contain a certain I element
    def find_indexes_with_element(self, element):
        indexes = []

        if element not in self.Items:
            return []

        element_idx = self.Items.index(element)
        n = len(self.Items)

        for combination in range(2 ** (n - 1)):

            binary_repr = format(combination, f'0{n - 1}b')

            new_binary_repr = binary_repr[:element_idx] + '1' + binary_repr[element_idx:]

            indexes.append(int(new_binary_repr, 2))

        return indexes

    # Returns all possible non-dominated N_c given an i in I_c and a w in W
    def n_c_generator(self, I_c, _w_max):

        widths = [item[0] for item in I_c]

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

    # returns kmin_{cj}, kmax_{cj}, lmin_{cjk}
    def optimised_stocksize_variables(self, index, _width):
        subset = self.powerset_at_index(index)
        _demands = [item[2] for item in subset]
        # lmin_{cjk}
        best_minimum_stock_size_length = {}

        for n_c in self.n_c_generator(subset, _width):
            # ceiling(a/b) = -floor(a/-b)
            # c_{ijn_i}
            columns_to_satisfy_demand = [-(_demands[i] // -n_c[i]) for i in range(len(_demands))]
            # cmax_{i}
            maximum_columns = [self.L_max // item[1] for item in subset]
            # kmin_{ijn_i}
            _minimum_stock_sizes = [-(columns_to_satisfy_demand[i] // -maximum_columns[i])
                                    for i in range(len(columns_to_satisfy_demand))]
            # cmin_{i}
            minimum_columns = [max(self.L_min, item[1]) // item[1] for item in subset]
            # kmax_{ijn_i}
            _maximum_stock_sizes = [-(columns_to_satisfy_demand[i] // -minimum_columns[i])
                                    for i in range(len(columns_to_satisfy_demand))]

            # kmin_cj
            _minimum_stock_size = max(_minimum_stock_sizes)
            # kmax_cj
            _maximum_stock_size = max(_maximum_stock_sizes)

            for k in range(_minimum_stock_size, _maximum_stock_size + 1):
                # chat_{ijn_ik}
                columns_per_panel = [-(column_to_satisfy_item_demand // -k)
                                    for column_to_satisfy_item_demand in columns_to_satisfy_demand]
                item_minimum_stock_size_length = [subset[i][1] * columns_per_panel[i] for i in range(len(subset))]
                maximum_item_minimum_stock_size_length = max(item_minimum_stock_size_length)
                if maximum_item_minimum_stock_size_length > self.L_max:
                    print(
                        f"Infeasible to meet the demand of I_{index} with {k} panels of width {_width} with row distribution of {n_c}")
                else:
                    if best_minimum_stock_size_length.get((index, self.available_widths.index(_width), k),
                                                        32767) > maximum_item_minimum_stock_size_length > self.L_min:
                        best_minimum_stock_size_length[(index, self.available_widths.index(_width), k)] = maximum_item_minimum_stock_size_length
                        best_minimum_stock_size = _minimum_stock_size
                        best_maximum_stock_size = _maximum_stock_size
                    elif maximum_item_minimum_stock_size_length < self.L_min:
                        best_minimum_stock_size_length[(index, self.available_widths.index(_width), k)] = self.L_min
                        best_minimum_stock_size = _minimum_stock_size
                        best_maximum_stock_size = _maximum_stock_size
        return best_minimum_stock_size, best_maximum_stock_size, best_minimum_stock_size_length

    # 4. Define subsets of I compatible with stock size j (C_j)
    # Although referring to the J matrix, it only depends on W values.
    # J[idx, :] --> C_j[idx]
    def C_j_generator(self):
        C_j = [[] for _ in range(len(self.available_widths))]
        for idx, width in enumerate(self.available_widths):
            for index in range(1, 2 ** len(self.Items) + 1):
                # Check that n_i_max is respected
                if bin(index).count('1') > self.n_i_max:
                    continue
                subset = self.powerset_at_index(index)
                # Check that One Group policy is respected
                if self.one_group and any(item != subset[0][1] for item in subset):
                    continue
                total_width = sum(item[0] for item in subset)
                if total_width <= width:
                    C_j[idx].append(index)
        return C_j




