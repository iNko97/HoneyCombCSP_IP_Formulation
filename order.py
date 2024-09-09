import numpy as np
from itertools import product, combinations
import pandas as pd

def is_dominated(new_combo, valid_combinations):
    for combination in valid_combinations:
        # Check if combo dominates new_combo
        if all(combination[i] >= new_combo[i] for i in range(len(new_combo))):
            return True
    return False

def remove_dominated(valid_combinations, new_combo):
    return [combination for combination in valid_combinations
            if not all(new_combo[i] >= combination[i] for i in range(len(combination)))]


# Returns all possible non-dominated N_c given an i in I_c and a w in W
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


class Order:
    def __init__(self, datasheet, scenario, order_number):
        self.order_sheet = pd.read_excel(datasheet,
                                         engine="odf",
                                         sheet_name=f'O{order_number}',
                                         dtype={'Item': None, 'Width': np.int16, 'Length': np.int16, 'Demand': np.int16}
                                         )
        self.param_sheet = pd.read_excel(datasheet,
                                         engine="odf",
                                         sheet_name=f'scenario_{scenario[0]}_{scenario[1]}'
                                         )
        print(self.param_sheet)

        # Factory Settings
        self.L_min = int(self.param_sheet['lmin'][0])   # Minimum panel length
        self.L_max = int(self.param_sheet['lmax'][0])  # Maximum panel length
        self.n_s_max = int(self.param_sheet['n_s_max'][0])  # Maximum number of stock sizes
        self.n_w_max = int(self.param_sheet['n_w_max'][0])  # Maximum number of widths
        self.n_i_max = int(self.param_sheet['n_i_max'][0])  # Maximum number of items per pattern
        self.one_group = bool(self.param_sheet['one_group'][0])  # Only one-groups are allowed
        self.available_widths = self.param_sheet['widths'].to_list()  # Set of w available stock widths
        print(self.available_widths)

        # Order information
        self.item_widths = self.order_sheet['Width']  # Individual item width
        self.item_lengths = self.order_sheet['Length']  # Individual item Length
        self.item_demands = self.order_sheet['Demand']  # Individual item Demand
        self.Items = [[int(self.item_widths[i]), int(self.item_lengths[i]), int(self.item_demands[i])]
                      for i in range(len(self.item_widths))]

        # Initial Pre-processing
        self.best_nc = {}
        print("Generating C")
        self.C = set()
        self.C_generator()
        print("Generating a_ic")
        self.a_ic = {}
        self.a_ic_generator()
        self.A_c = {}

    # Returns the items subset I_c corresponding to a specified index.
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

    # returns kmin_{cj}, kmax_{cj}, lmin_{cjk}
    def optimised_stocksize_variables(self, index, _width):
        subset = self.powerset_at_index(index)
        _demands = [item[2] for item in subset]
        lmin_cjk = {}
        K_cj = set()

        for n_c in n_c_generator(subset, _width):
            # ceiling(a/b) = -floor(a/-b)
            c_ijni = [-(_demands[i] // -n_c[i]) for i in range(len(_demands))]

            cmax_i = [self.L_max // item[1] for item in subset]
            cmin_i = [max(self.L_min, item[1]) // item[1] for item in subset]

            kmin_ijni = [-(c_ijni[i] // -cmax_i[i]) for i in range(len(c_ijni))]
            kmax_ijni = [-(c_ijni[i] // -cmin_i[i]) for i in range(len(c_ijni))]

            kmin_cj = max(kmin_ijni)
            kmax_cj = max(kmax_ijni)

            K_cj.update([k for k in range(kmin_cj, kmax_cj+1)])

            for k in range(kmin_cj, kmax_cj + 1):
                chat_ijnik = [-(c // -k) for c in c_ijni]
                l_icjnik = [subset[i][1] * chat_ijnik[i] for i in range(len(subset))]
                lmin_Nccjk = max(max(l_icjnik), self.L_min)

                if lmin_Nccjk > self.L_max:
                    print(
                        f"Infeasible to meet the demand of I_{index} with {k} panels "
                        f"of width {_width} with row distribution of {n_c}"
                    )
                else:
                    if lmin_cjk.get((index, self.available_widths.index(_width), k), lmin_Nccjk+1) > lmin_Nccjk:
                        lmin_cjk[(index, self.available_widths.index(_width), k)] = lmin_Nccjk
                        self.best_nc[(index, self.available_widths.index(_width), k)] = n_c

        return K_cj, lmin_cjk

    # 4. Define subsets of I compatible with stock size j (C_j)
    # Although referring to the J matrix, it only depends on W values.
    # J[idx, :] --> C_j[idx]
    def C_j_generator(self):
        C_j = [[] for _ in range(len(self.available_widths))]
        for idx, width in enumerate(self.available_widths):
            for index in self.C:
                subset = self.powerset_at_index(index)
                # Check that One Group policy is respected
                if self.one_group and any(item[1] != subset[0][1] for item in subset):
                    continue

                total_width = sum(item[0] for item in subset)

                if total_width <= width:
                    C_j[idx].append(index)
        return C_j

    def C_generator(self):
        for max_items in range(1, self.n_i_max+1):
            n = 2**(len(self.Items))-1
            max_bits = n.bit_length()
            for num_bits in range(max_items, max_bits + 1):
                for ones_positions in combinations(range(num_bits), max_items):
                    num = sum(1 << pos for pos in ones_positions)
                    if num > n:
                        break
                    self.C.add(num)
        self.C = sorted(self.C)

    # Checks whether item i is in index c
    def a_ic_generator(self):
        n = len(self.Items)
        for c in self.C:
            for i in range(n):  # 0 n-1
                if bool(c & (1 << i)):
                    self.a_ic[(n-i-1, c)] = 1
                else:
                    self.a_ic[(n-i-1, c)] = 0
