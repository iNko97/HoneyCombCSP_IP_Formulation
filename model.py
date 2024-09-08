import gurobipy as gp
from gurobipy import GRB
from order import Order
import csv, os

def optimise(order_number, scenario_id, _n_s_max):
    # Program Parameters
    path = "./Data/Input_data.ods"
    order_filename = f"./Output/order_{order_number}_{scenario_id}_{_n_s_max}.csv"
    results_filename = f"./Output/results.csv"
    write_header = not os.path.exists(results_filename)
    item_allocations = []
    A_c_lower_bound = 0

    # Gurobi Parameters
    model = gp.Model("2D_Cutting_Stock")
    model.setParam(GRB.Param.MIPFocus, 3)
    model.setParam(GRB.Param.PreDual, -1)
    model.setParam(GRB.Param.MIPGap, 0)
    model.setParam(GRB.Param.TimeLimit, 1800)
    model.setParam(GRB.Param.Heuristics, 1)
    model.setParam(GRB.Param.Presolve, 2)
    model.setParam(GRB.Param.Cuts, 3)

    # Support variables

    # Factory settings
    order = Order(path, (scenario_id, _n_s_max), order_number)

    L_min = order.L_min  # Minimum panel length
    L_max = order.L_max  # Maximum panel length
    n_s_max = order.n_s_max  # Maximum number of stock sizes
    n_w_max = order.n_w_max  # Maximum number of widths
    # n_i_max = order.n_i_max  # Maximum number of items per pattern
    # one_group = order.one_group  # Only one-groups are allowed
    W = order.available_widths  # Set of w available stock widths
    I = order.Items  # I: item types with their Width, Length, and Demand

    # SETS AND PARAMETERS

    # Potential stocks J of dimensions |W| * n_s_max
    # This is essentially a matrix of all x_j
    # Each row is indexed by available widths.
    # For a certain width w, if the index idx = W.index(w),then the set J^w' = J[idx, :]
    J = (len(W), n_s_max)
    a_ic = order.a_ic

    # Subsets of I compatible with stock size j (C_j)
    # Although referring to the J matrix, it only depends on W values.
    # given a j at index J[idx, 0], ..., J[idx, n] --> C_j[idx] forall n
    print("Generating C_j")
    C_j = order.C_j_generator()

    print("Generating K_cj and lmin_cjk")
    # dictionary with triplet (c, idx_j, k) where c in I_c, idx_j in J.shape[0], k stock size
    lmin_cjk = {}
    # dictionary with tuple (c, idx_j) : list(range(kmin_cj, ..., kmax_cj))
    K_cj = {}
    #  dictionary with (idx) : Delta_j, that is the minimum difference between two lmin_cjk
    Delta_j = {}

    for idx in range(J[0]):
        for c in C_j[idx]:
            (_K_cj, _lmin_cjk) = order.optimised_stocksize_variables(c, W[idx])
            lmin_cjk.update(_lmin_cjk)
            K_cj[(c, idx)] = _K_cj

    # PRE-PROCESSING AND MODEL DIMENSIONS REDUCTION

    # Proposition 1: Generalised not just for k-1 but for the biggest k_prime smaller than k
    # Instead of taking k-1, we take the maximum k_prime lower than k. This allows for tighter lower bounds.
    print("Applying Proposition 1.")
    for idx in range(J[0]):
        for c in C_j[idx]:
            k_to_be_removed = []
            for k in K_cj[c, idx]:
                if k == min(K_cj[c, idx]):
                    continue
                k_prev = max((x for x in K_cj[c, idx] if x < k))
                if lmin_cjk[c, idx, k_prev] == lmin_cjk[c, idx, k]:
                    k_to_be_removed.append(k)
            for k in k_to_be_removed:
                K_cj[c, idx].remove(k)

    # Proposition 2: Substituting I_c with |I_c| cutting pattern with lower lower bounds
    print("Applying Proposition 2.")
    for idx in range(J[0]):
        for c in C_j[idx]:
            k_to_be_removed = []
            c_primes = sorted([1 << i for i in range(c.bit_length()) if c & (1 << i)], reverse=True)
            if len(c_primes) == 1:
                continue
            for k in K_cj[c, idx]:
                k_i = {}
                for c_prime in c_primes:
                    for k_prime in K_cj[c_prime, idx]:
                        if lmin_cjk[c_prime, idx, k_prime] <= lmin_cjk[c, idx, k]:
                            k_i[c_prime] = min(k_i.get(c_prime, k_prime), k_prime)
                if len(k_i) == len(c_primes) and sum(k_i.values()) <= k:
                    k_to_be_removed.append(k)
            for k in k_to_be_removed:
                K_cj[c, idx].remove(k)

    # \Delta_j, the minimum length difference between two different values of lmin_cjk for a stock j
    # Instead of taking different c_1 and c_2 as suggested, we take c_1 and c_2 with no items in common, as an item may
    # not be assigned to multiple sock sizes. This allows for larger steps and diminution of the search space.
    print("Generating Delta_j.")
    for idx in range(J[0]):
        best_Delta = L_max
        for c_1 in C_j[idx]:
            for c_2 in C_j[idx]:
                if c_1 & c_2 > 0.5:
                    continue
                K_cj_1 = K_cj[(c_1, idx)]
                K_cj_2 = K_cj[(c_2, idx)]
                for k_1 in K_cj_1:
                    for k_2 in K_cj_2:
                        diff = abs(lmin_cjk[(c_1, idx, k_1)] - lmin_cjk[(c_2, idx, k_2)])
                        if 0 < diff < best_Delta:
                            best_Delta = diff
        Delta_j[idx] = best_Delta

    # A_c lower bound for the area of material used
    print("Calculating A_c lower bounds.")
    A_c = {}
    for (c, idx, k), _lmin_cjk in lmin_cjk.items():
        current_product = k * W[idx] * _lmin_cjk
        A_c[c] = min(A_c.get(c, current_product), current_product)
    order.A_c = A_c

    n_c_asterisk = order.best_nc

    print("Initialising model.")
    # DECISION VARIABLES

    # \alpha_{cj} 1 if the subset of item types I_c is assigned to stock size j, 0 otherwise for j \in J, c \in C
    alpha_cj = model.addVars(
        [(c, idx, n)
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]
         ],
        vtype=GRB.BINARY,
        name="alpha_cj"
    )

    # \beta_j 1 if stock size j is selected, 0 otherwise for j \in J
    beta_j = model.addVars(
        [(idx, n)
         for idx in range(J[0])
         for n in range(J[1])
         ],
        vtype=GRB.BINARY,
        name="beta_j"
    )

    # \gamma_{cjk} is 1 if a total of k panels of stock size j are produced to satisfy the demand for item types in
    # subset I_c, and is 0 otherwise, for j \in J, c \in C_j, k \in K_{cj}
    gamma_cjk = model.addVars(
        [(c, idx, n, k)
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]
         for k in K_cj[(c, idx)]
         ],
        vtype=GRB.BINARY,
        name="gamma_cjk"
    )

    # \delta_{w} 1 if at least one stock size in J^{w} is selected, 0 otherwise for w_prime \in W
    delta_w = model.addVars(
        [w
         for w in range(len(W))
         ],
        vtype=GRB.BINARY,
        name="delta_w"
    )

    # x_j length of stock size j, for j \in J
    x_j = model.addVars(
        [(idx, n)
         for idx in range(J[0])
         for n in range(J[1])
         ],
        vtype=GRB.CONTINUOUS,
        ub=L_max,
        name="x_j"
    )

    # y_{cjk}
    # Auxiliary continuous variable needed to linearise the product x_j*\gamma_{cjk} for j \in J, c \in C, k \in K_{cj}.
    # Note that y_{cjk} = x_j if \gamma_{cjk} = 1, otherwise 0
    y_cjk = model.addVars(
        [(c, idx, n, k)
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]
         for k in K_cj[(c, idx)]],
        vtype=GRB.CONTINUOUS,
        lb=0,
        ub=L_max,
        name="y_cjk"
    )

    # OBJECTIVE FUNCTION
    # 2. Minimise the total area of the material used
    model.setObjective(
        gp.quicksum(
            W[idx] * k * y_cjk[c, idx, n, k]
            for idx in range(J[0])
            for n in range(J[1])
            for c in C_j[idx]
            for k in K_cj[(c, idx)]
        ),
        GRB.MINIMIZE
    )

    # CONSTRAINTS

    # 3. Ensure that each item is allocated to a stock size
    # noinspection PyTypeChecker
    model.addConstrs(
        (gp.quicksum(alpha_cj[c, idx, n] * a_ic[(i, c)]
                     for idx in range(J[0])
                     for n in range(J[1])
                     for c in C_j[idx]) == 1
         for i in range(len(I))),
        name='link_alpha_a'
    )

    # 4. Ensure that a stock size is used iff at least 1 item type is assigned to it
    model.addConstrs(
        (alpha_cj[c, idx, n] <= beta_j[idx, n]
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]),
        name="link_alpha_beta"
    )

    # 4b.
    model.addConstrs(
        (gp.quicksum(alpha_cj[c, idx, n]
                     for c in C_j[idx]) >= beta_j[idx, n]
         for idx in range(J[0])
         for n in range(J[1])),
        name="link_beta_alpha"
    )

    # 5. Ensure the limit on stock sizes is not broken
    model.addConstr(
        gp.quicksum(
            beta_j[idx, n]
            for idx in range(J[0])
            for n in range(J[1])
        ) <= n_s_max,
        name="limit_stock_sizes"
    )

    # 6. Ensure that if a subset I_c is assigned to j, then a k is assigned.
    model.addConstrs(
        (gp.quicksum(gamma_cjk[c, idx, n, k] for k in K_cj[(c, idx)]) == alpha_cj[c, idx, n]
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]),
        name="link_gamma_alpha"
    )

    # 7. Ensure y_cjk and x_j constraints
    model.addConstrs(
        (gp.quicksum(lmin_cjk[(c, idx, k)] * gamma_cjk[c, idx, n, k]
                     for k in K_cj[(c, idx)]) <= x_j[idx, n]
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]),
        name="length_bound_7"
    )

    # 8. Ensure that the y_cjk respects its length st.
    model.addConstrs(
        (y_cjk[c, idx, n, k] <= x_j[idx, n]
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]
         for k in K_cj[(c, idx)]),
        name="link_8_y_x"
    )

    # 9. Ensure that y_cjk doesn't exceed maximum length
    model.addConstrs(
        (y_cjk[c, idx, n, k] <= L_max * gamma_cjk[c, idx, n, k]
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]
         for k in K_cj[(c, idx)]),
        name="link_9_y_x"
    )

    # 10. Ensure that x_j doesn't exceed maximum length if k panels are being produced.
    model.addConstrs(
        (x_j[idx, n] - y_cjk[c, idx, n, k] <= L_max * (1 - gamma_cjk[c, idx, n, k])
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]
         for k in K_cj[(c, idx)]),
        name="link_10_y_x"
    )

    # 11. Ensure x_j doesn't exceed minimum length for all selected stock sizes j, otherwise 0.
    model.addConstrs(
        (L_min * beta_j[idx, n] <= x_j[idx, n]
         for idx in range(J[0])
         for n in range(J[1])),
        name="stock_length_min"
    )
    # 11b.
    model.addConstrs(
        (x_j[idx, n] <= L_max * beta_j[idx, n]
         for idx in range(J[0])
         for n in range(J[1])),
        name="stock_length_max"
    )

    # 12. Ensure a stock size width is used iff used at least one time.
    model.addConstrs(
        (gp.quicksum(beta_j[idx, n] for n in range(J[1])) >= delta_w[idx]
         for idx in range(J[0])),
        name="link_delta_beta"
    )

    # 12b.
    model.addConstrs(
        (beta_j[idx, n] <= delta_w[idx]
         for idx in range(J[0])
         for n in range(J[1])),
        name="link_beta_delta"
    )

    # 13. Limit the number of stock widths
    model.addConstr(
        gp.quicksum(delta_w[idx] for idx in range(J[0])) <= n_w_max,
        name="limit_stock_widths"
    )

    # 21. Symmetry breaking constraint for \beta_j
    model.addConstrs(
        (beta_j[idx, n + 1] <= beta_j[idx, n]
         for idx in range(J[0])
         for n in range(J[1] - 1)),
        name="beta_symmetry"
    )

    # 22. Symmetry breaking constraints for \x_j
    model.addConstrs(
        (x_j[idx, n] - x_j[idx, n + 1] >= Delta_j[idx] * beta_j[idx, n + 1]
         for idx in range(J[0])
         for n in range(J[1] - 1)),
        name="x_symmetry"
    )

    # Custom. Cj sum is always the same
    model.addConstr(
        gp.quicksum(c * alpha_cj[c, idx, n]
                    for idx in range(J[0])
                    for n in range(J[1])
                    for c in C_j[idx]) == 2 ** (len(I)) - 1,
        name="link_max_sum"
    )

    # Optimize the model
    model.optimize()

    # Extract the solution
    if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
        if model.status == GRB.OPTIMAL:
            print("Optimal solution found:")
        if model.status == GRB.TIME_LIMIT:
            print("Time limit reached. Best solution found:")
        for idx in range(J[0]):
            for n in range(J[1]):
                if beta_j[idx, n].x > 0.5:
                    print(f"Stock size {idx}, {n}: width = {W[idx]}, length = {x_j[idx, n].x}")
                    for c in C_j[idx]:
                        if alpha_cj[c, idx, n].x > 0.5:
                            binary_rep = f"{c:0{len(I)}b}"
                            indexes = [i + 1 for i, bit in enumerate(binary_rep) if bit == '1']
                            A_c_lower_bound += order.A_c[c]
                            for k in K_cj[(c, idx)]:
                                if gamma_cjk[c, idx, n, k].x > 0.5:
                                    print(f"    {k} panels of items {indexes}")
                                    print(f"        with respectively {n_c_asterisk[(c, idx, k)]} columns.")
                                    results.append({
                                        "Stock_Length": W[idx],
                                        "Stock_Width": x_j[idx, n].x,
                                        "Panels": k,
                                        "Items": indexes,
                                        "Columns": n_c_asterisk[(c, idx, k)]
                                    })
        print(f"The lower bound for this solution instance was {round(A_c_lower_bound/1000000)}")

        with open(order_filename, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=["Stock_Length", "Stock_Width", "Panels", "Items", "Columns"])
            writer.writeheader()
            for _item in item_allocations:
                writer.writerow(_item)

        with open(results_filename, mode='a', newline='') as file:
            if write_header:
                writer = csv.DictWriter(file,
                                        fieldnames=["Order", "Scenario", "n_s_max", "Solution", "Gap", "Time", "A_c"])
                writer.writeheader()
            result = {"Order": order_number,
                      "Scenario": scenario_id,
                      "n_s_max": n_s_max,
                      "Solution": round(model.getObjective()/1000000),
                      "Gap": model.MIPGap,
                      "Time": model.Runtime,
                      "A_c": round(A_c_lower_bound/1000000)}
            writer.writerow(result)


    else:
        print("No solution found.")
    return model, order
