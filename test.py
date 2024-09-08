import csv
from model import optimise

# Define your ranges
order_numbers = [7]
scenario_ids = [1, 2, 3]
n_s_max_ids = [1, 2]

results = []

for order_number in order_numbers:
    for scenario_id in scenario_ids:
        for n_s_max_id in n_s_max_ids:

            model = optimise(order_number, scenario_id, n_s_max_id)

            for idx in range(J[0]):
                for n in range(J[1]):
                    if beta_j[idx, n].x > 0.5:
                        stock_size = f"Stock size {idx}, {n}: width = {W[idx]}, length = {x_j[idx, n].x}"

                        for c in C_j[idx]:
                            if alpha_cj[c, idx, n].x > 0.5:
                                binary_rep = f"{c:0{len(I)}b}"
                                indexes = [i + 1 for i, bit in enumerate(binary_rep) if bit == '1']

                                for k in K_cj[(c, idx)]:
                                    if gamma_cjk[c, idx, n, k].x > 0.5:
                                        panels = f"{k} panels of items {indexes}"
                                        columns = f"with respectively {n_c_asterisk[(c, idx, k)]} columns."

                                        # Append the result as a dictionary
                                        results.append({
                                            "Stock_Size": stock_size,
                                            "Panels": panels,
                                            "Columns": columns
                                        })

# Define the CSV file name
csv_filename = "optimisation_results.csv"

# Write the results to a CSV file
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["Stock_Size", "Panels", "Columns"])
    writer.writeheader()  # Write the header
    for result in results:
        writer.writerow(result)  # Write each row

print(f"Results saved to {csv_filename}")