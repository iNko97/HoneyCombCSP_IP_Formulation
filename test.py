# Change import between chained_LB_model, chained_model or model
# Model is the default model
# Chained model sets the initial solution of the next model as the optimal solution of the previous
# Chained LB also sets the lower bound of the optimal solution as A_c
from chained_model import optimise
import os

# Available orders: 7, 8, 10, 17, 20, 22, 23, 24
order_numbers = [10, 17, 20, 22, 23, 24]

# Available scenarios: 2-item 1, 2, 3 and respective 3-items 4, 5, 6
scenario_ids = [1, 2, 3, 6]

# List of directories to create
dirs = ["Gurobi", "Gurobi/chained_LB_model", "Gurobi/chained_model", "Output", "Output/chained_LB_model", "Output/chained_model", "Output/model"]

for directory in dirs:
    os.makedirs(directory, exist_ok=True)

for order_number in order_numbers:
    for scenario_id in scenario_ids:
        if scenario_id == 1:
            n_s_max = [1, 2]
        else:
            n_s_max = [2]
        for n_s_max_id in n_s_max:
            print(f"Optimising order {order_number} with scenario {scenario_id} and n_s_max_id {n_s_max_id}")
            _m, _o = optimise(order_number, scenario_id, n_s_max_id)
