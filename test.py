import csv
from model import optimise

# Define your ranges
order_numbers = [7]
scenario_ids = [1, 2, 3]
n_s_max_ids = [1, 2]

for order_number in order_numbers:
    for scenario_id in scenario_ids:
        for n_s_max_id in n_s_max_ids:
            print(f"Optimising order {order_number} with scenario {scenario_id} and n_s_max_id {n_s_max_id}")
            model, order = optimise(order_number, scenario_id, n_s_max_id)
