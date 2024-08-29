import numpy as np
from itertools import product
import pandas as pd

# Factory settings
#TODO: IMPORT FROM CSV

class Order:
    L_min = 700  # Minimum panel length
    L_max = 3100  # Maximum panel length
    n_s_max = 2  # Maximum number of stock sizes
    n_w_max = 2  # Maximum number of widths
    W = [1200, 1400, 1550, 1600]  # Set of w available stock widths
    I=[]
    def __init__(self, datasheet):
        number = input("Scegli il numero dell'ordine (da 1 a 8)")
        self.sheet = pd.read_excel(datasheet,number)
        self.Width = self.sheet['Width']
        self.Length = self.sheet['Length']
        self.Demand = self.sheet['Demand']
        self.I = [self.Width, self.Length,self.Demand]

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
    def n_c_generator(I_c, _w_max):

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

    def is_dominated(new_combo, valid_combinations):
        for combination in valid_combinations:
            # Check if combo dominates new_combo
            if all(combination[i] >= new_combo[i] for i in range(len(new_combo))):
                return True
        return False

    def remove_dominated(valid_combinations, new_combo):
        return [combination for combination in valid_combinations
                if not all(new_combo[i] >= combination[i] for i in range(len(combination)))]

    # returns kmin_{cj}, kmax_{cj}, lmin_{cjk}
    def optimised_stocksize_variables(index, _width):
        subset = powerset_at_index(index)
        _demands = [item[2] for item in subset]
        # lmin_{cjk}
        best_minimum_stock_size_length = {}

        for n_c in n_c_generator(subset, _width):
            # ceiling(a/b) = -floor(a/-b)
            # c_{ijn_i}
            columns_to_satisfy_demand = [-(_demands[i] // -n_c[i]) for i in range(len(_demands))]
            # cmax_{i}
            maximum_columns = [L_max // item[1] for item in subset]
            # kmin_{ijn_i}
            _minimum_stock_sizes = [-(columns_to_satisfy_demand[i] // -maximum_columns[i])
                                    for i in range(len(columns_to_satisfy_demand))]
            # cmin_{i}
            minimum_columns = [max(L_min, item[1]) // item[1] for item in subset]
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
                if maximum_item_minimum_stock_size_length > L_max:
                    print(
                        f"Infeasible to meet the demand of I_{index} with {k} panels of width {_width} with row distribution of {n_c}")
                else:
                    if best_minimum_stock_size_length.get((index, W.index(_width), k),
                                                        32767) > maximum_item_minimum_stock_size_length > L_min:
                        best_minimum_stock_size_length[(index, W.index(_width), k)] = maximum_item_minimum_stock_size_length
                        best_minimum_stock_size = _minimum_stock_size
                        best_maximum_stock_size = _maximum_stock_size
                    elif maximum_item_minimum_stock_size_length < L_min:
                        best_minimum_stock_size_length[(index, W.index(_width), k)] = L_min
                        best_minimum_stock_size = _minimum_stock_size
                        best_maximum_stock_size = _maximum_stock_size
        return best_minimum_stock_size, best_maximum_stock_size, best_minimum_stock_size_length


    # 3. Define potential stocks (J) of dimensions W * n_s_max
    # each row indexes are ordered by available widths.
    # J = np.zeros((len(W), n_s_max))

    # 4. Define subsets of I compatible with stock size j (C_j)
    # Although referring to the J matrix, it only depends on W values.
    # J[idx, :] --> C_j[idx]
    def C_j_generator():
        C_j = [[] for _ in range(len(W))]
        for idx, width in enumerate(W):
            for index in range(1, 2**len(I)+1):
                subset = powerset_at_index(index)
                total_width = sum(item[0] for item in subset)

                if total_width <= width:
                    C_j[idx].append(index)
        return C_j

    #Checks whether item i is in index c
    def a_ic_generator(i, c):
        return 1 if bool(c & (1 << i)) else 0


