from model import optimise

# Gli scenari descritti sono gli stessi del paper.
# Si possono accedere alle soluzioni calcolate nei CSV nella cartella di output:
#   results.csv fornisce le soluzioni equivalenti alla tabella 5 del paper.
#   I vari files order_a_b_c sono le effettive soluzioni calcolate,
#       a e' il numoro dell'ordine
#       b e' il numero dello scenario
#       c e' il valore massimo di stock definibili in quello scenario

order_numbers = [7]
scenario_ids = [1]
n_s_max_ids = [1]

for order_number in order_numbers:
    for scenario_id in scenario_ids:
        for n_s_max_id in n_s_max_ids:
            if scenario_id == 3 and n_s_max_id == 1:
                continue
            print(f"Optimising order {order_number} with scenario {scenario_id} and n_s_max_id {n_s_max_id}")
            _m, _o = optimise(order_number, scenario_id, n_s_max_id)
