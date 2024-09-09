import csv
import pandas as pd
import numpy as np

# Order to be verified for feasibility.
order_number = 24
scenario_number = 3
n_s_max = 2


def read_solution_from_csv(file_path):
    solution = []
    with open(file_path, mode='r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)
        for row in csv_reader:
            stock_length = int(row[0])
            stock_width = int(row[1])
            panels = int(row[2])
            items_list = eval(row[3])
            columns = eval(row[4])
            solution.append((stock_length, stock_width, panels, items_list, columns))
    return solution


def check_constraints(items, solution):
    error = False
    item_dict = {i+1: {'width': item[0], 'length': item[1], 'demand': item[2]} for i, item in
                 enumerate(items)}

    for stock_width, stock_length, panels, items_list, columns in solution:
        total_width = 0
        for item_index, column_count in zip(items_list, columns):
            item_width = item_dict[item_index]['width']
            item_length = item_dict[item_index]['length']

            total_width += column_count * item_width

            items_produced = (stock_length // item_length) * column_count * panels
            if items_produced < item_dict[item_index]['demand']:
                print(
                    f"ERROR Demand constraint violated for item {item_index}: {items_produced} < Demand {item_dict[item_index]['demand']}")
                error = True
            else:
                print(f"    Item {item_index} produced {items_produced} and Satisfied demand of {item_dict[item_index]['demand']}")

        if total_width > stock_width:
            print(f"ERROR Stock width exceeded for {items_list}, {columns}")
            error = True
        else:
            print(f"    Width OK for {items_list}, {columns}, for {total_width} <= {stock_width}")
    return error


order_sheet = pd.read_excel(f"./Data/Input_data.ods",
                            engine="odf",
                            sheet_name=f'O{order_number}',
                            dtype={'Item': None, 'Width': np.int16, 'Length': np.int16, 'Demand': np.int16}
                            )

# Items information
item_widths = order_sheet['Width']
item_lengths = order_sheet['Length']
item_demands = order_sheet['Demand']
items = [[int(item_widths[i]), int(item_lengths[i]), int(item_demands[i])] for i in range(len(item_widths))]

csv_file_path = f'./Output/order_{order_number}_{scenario_number}_{n_s_max}.csv'

solution = read_solution_from_csv(csv_file_path)

# Check constraints
print(f"Errors? {check_constraints(items, solution)}")


