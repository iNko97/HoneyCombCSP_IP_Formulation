from model import optimise

order_numbers = [8, 10, 24]
scenario_ids = [1, 2, 3]
n_s_max_ids = [1, 2]

for order_number in order_numbers:
    for scenario_id in scenario_ids:
        for n_s_max_id in n_s_max_ids:
            if scenario_id == 3 and n_s_max_id == 1:
                continue
            print(f"Optimising order {order_number} with scenario {scenario_id} and n_s_max_id {n_s_max_id}")
            _m, _o = optimise(order_number, scenario_id, n_s_max_id)