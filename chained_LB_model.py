import csv
import os
import gurobipy as gp
from gurobipy import GRB

from random import shuffle
from order import Order
from collections import Counter


def optimise(order_number, scenario_id, _n_s_max):
    # Program Parameters
    path = "./Data/Input_data.ods"
    order_filename = f"./Output/chained_LB_model/order_{order_number}_{scenario_id}_{_n_s_max}.csv"
    results_filename = f"./Output/chained_LB_model/results.csv"
    write_header = not os.path.exists(results_filename)
    item_allocations = []

    # Gurobi Parameters
    model = gp.Model("2D_Cutting_Stock")
    model.setParam(GRB.Param.TimeLimit, 1800)
    model.setParam(GRB.Param.MIPGap, 0.000001)
    model.setParam(GRB.Param.MIPFocus, 2)
    model.setParam(GRB.Param.Presolve, 2)
    model.setParam(GRB.Param.Cuts, 3)
    model.setParam(GRB.Param.Method, 2)
    model.setParam('NumericFocus', 2)

    # Factory settings
    order = Order(path, (scenario_id, _n_s_max), order_number)
    n_s_max = order.n_s_max  # Maximum number of stock sizes
    n_w_max = order.n_w_max  # Maximum number of widths
    W = order.available_widths  # Set of w available stock widths
    I = order.Items  # I: item types with their Width, Length, and Demand

    # SETS AND PARAMETERS

    # Potential stocks J of dimensions |W| * n_s_max
    # This is essentially a matrix of all x_j
    # Each row is indexed by available widths.
    # For a certain width w, if the index idx = W.index(w),then the set J^w' = J[idx, :]
    J = (len(W), n_s_max)

    # Whether an item i belongs to I_c
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
            k_min = min(K_cj[c, idx])
            for k in K_cj[c, idx]:
                if k == k_min:
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

    # A_c lower bound for the area of material used
    print("Calculating A_c lower bounds.")
    A_c = {}
    lower_bound = 0

    ac_model = gp.Model("Unbound_Cutting_Stock")
    for (c, idx, k), _lmin_cjk in lmin_cjk.items():
        current_product = k * W[idx] * _lmin_cjk
        A_c[c] = min(A_c.get(c, current_product), current_product)
    order.A_c = A_c

    print("Initialising Ac model.")

    # \alpha_{c} 1 if the subset of item types I_c is assigned to stock size j, 0 otherwise for j \in J, c \in C
    ac_alpha_c = ac_model.addVars(
        [c
         for c in C_j[idx]
         ],
        vtype=GRB.BINARY,
        name="ac_alpha_c"
    )

    #Objective for A_c model
    ac_model.setObjective(
        gp.quicksum(
            A_c[c] * ac_alpha_c[c]
            for c in C_j[idx]
        ),
        GRB.MINIMIZE
    )

    #s.t. for A_c model
    ac_model.addConstrs(
        (gp.quicksum(ac_alpha_c[c] * a_ic[(i, c)]
                     for c in C_j[idx]) == 1
         for i in range(len(I))),
        name='ac_link_alpha_a'
    )

    ac_model.optimize()

    # Extract the solution
    if ac_model.status == GRB.OPTIMAL:
        print(f"A_c Done! {ac_model.ObjVal}")
        lower_bound = ac_model.ObjVal




    # Get column distribution data for solution parsing.
    n_c_asterisk = order.best_nc

    # L_j matrix of available lengths.
    print("Generating L_j")
    stock_lengths = {}
    for idx in range(J[0]):
        stock_lengths[idx] = {}
        lista = sorted(list(Counter(lmin_cjk[c, idx, k] for c in C_j[idx] for k in K_cj[(c, idx)])))
        for i, length in enumerate(lista):
            stock_lengths[idx][i] = length
        print(stock_lengths[idx], len(stock_lengths[idx]))

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

    # \xi_jl index of l in stock_lengths[idx] length of stock size j, for j \in J
    xi_jl = model.addVars(
        [(idx, n, l)
         for idx in range(J[0])
         for n in range(J[1])
         for l in range(len(stock_lengths[idx]))
         ],
        vtype=GRB.BINARY,
        name="xi_jl"
    )

    # Artificial objective function z
    z = model.addVar(name="z", lb=lower_bound)

    # OBJECTIVE FUNCTION
    # 2. Minimise the total area of the material used
    model.setObjective(z, GRB.MINIMIZE)

    # CONSTRAINTS
    #2a. z equal to the original objective function value
    model.addConstr(
        z == gp.quicksum(
            W[idx] *
            k * gamma_cjk[c, idx, n, k] *
            gp.quicksum(xi_jl[idx, n, l] * stock_lengths[idx][l] for l in range(len(stock_lengths[idx])))
            for idx in range(J[0])
            for n in range(J[1])
            for c in C_j[idx]
            for k in K_cj[(c, idx)]
        ),
        name="z_constraint"
    )

    # 3. Ensure that each item is allocated to a stock size
    model.addConstrs(
        (gp.quicksum(alpha_cj[c, idx, n] * a_ic[(i, c)]
                     for idx in range(J[0])
                     for n in range(J[1])
                     for c in C_j[idx]) == 1
         for i in range(len(I))),
        name='link_alpha_a'
    )

    # New. Ensure that each stock size is allocated a length iff it's allocated.
    model.addConstrs(
        (gp.quicksum(xi_jl[idx, n, l] for l in range(len(stock_lengths[idx]))) == beta_j[idx, n]
         for idx in range(J[0])
         for n in range(J[1])),
        name='link_xi_jl'
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
        (gp.quicksum(alpha_cj[c, idx, n] for c in C_j[idx]) >= beta_j[idx, n]
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

    # 7. Ensure \gamma_cjk and xi_jl constraints
    model.addConstrs(
        (gp.quicksum(lmin_cjk[(c, idx, k)] * gamma_cjk[c, idx, n, k] for k in K_cj[(c, idx)]) <=
         gp.quicksum(xi_jl[idx, n, l] * stock_lengths[idx][l] for l in range(len(stock_lengths[idx])))
         for idx in range(J[0])
         for n in range(J[1])
         for c in C_j[idx]),
        name="length_bound_7"
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

    # 22. Symmetry breaking constraints for \xi_jl, we choose bigger indexes first
    model.addConstrs(
        (gp.quicksum((l+1) * xi_jl[idx, n, l] for l in range(len(stock_lengths[idx]))) -
         gp.quicksum((l+1) * xi_jl[idx, n+1, l] for l in range(len(stock_lengths[idx])))
         >= beta_j[idx, n+1]
         for idx in range(J[0])
         for n in range(J[1] - 1)),
        name="xi_symmetry"
    )

    # 23. Symmetry breaking st.
    # model.addConstrs(
    #     (gp.quicksum(xi_jl[idx, n, l] for n in range(J[1])) <= 1
    #      for idx in range(J[0])
    #      for l in range(len(stock_lengths[idx]))),
    #     name="xi_symmetry"
    # )

    model.update()

    # Load old solution
    if _n_s_max == 2:
        if scenario_id == 1:
            load_scenario_id = 1
            load_ns_max = 1
        elif scenario_id == 6:
            load_scenario_id = 3
            load_ns_max = 2
        else:
            load_scenario_id = scenario_id - 1
            load_ns_max = 2
        model.read(f"./Gurobi/chained_LB_model/solution_{order_number}_{load_scenario_id}_{load_ns_max}.mst")

    # Optimize the model
    model.optimize()

    # Extract the solution
    if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT or model.status == GRB.INTERRUPTED:
        if model.status == GRB.OPTIMAL:
            print("Optimal solution found:")
        if model.status == GRB.TIME_LIMIT or model.status == GRB.INTERRUPTED:
            print("Time limit reached or interrupted. Best solution found:")

        model.write(f"./Gurobi/chained_LB_model/solution_{order_number}_{scenario_id}_{_n_s_max}.mst")

        for idx in range(J[0]):
            for n in range(J[1]):
                if beta_j[idx, n].x > 0.5:
                    for l in range(len(stock_lengths[idx])):
                        if xi_jl[idx, n, l].x > 0.5:
                            print(f"Stock size {idx}, {n}: width = {W[idx]}, length = {stock_lengths[idx][l]}")
                            for c in C_j[idx]:
                                if alpha_cj[c, idx, n].x > 0.5:
                                    binary_rep = f"{c:0{len(I)}b}"
                                    indexes = [i + 1 for i, bit in enumerate(binary_rep) if bit == '1']
                                    for k in K_cj[(c, idx)]:
                                        if gamma_cjk[c, idx, n, k].x > 0.5:
                                            print(f"    {k} panels of items {indexes}")
                                            print(f"        with respectively {n_c_asterisk[(c, idx, k)]} columns.")
                                            item_allocations.append({
                                                "Stock_Length": W[idx],
                                                "Stock_Width": stock_lengths[idx][l],
                                                "Panels": k,
                                                "Items": indexes,
                                                "Columns": n_c_asterisk[(c, idx, k)]
                                            })
        print(f"The lower bound for this solution instance was {round(lower_bound / 1000000)}")

        with open(order_filename, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=["Stock_Length", "Stock_Width", "Panels", "Items", "Columns"])
            writer.writeheader()
            for _item in item_allocations:
                writer.writerow(_item)

        with open(results_filename, mode='a', newline='') as file:
            writer = csv.DictWriter(file,
                                    fieldnames=["Order", "Scenario", "n_s_max", "Solution", "Gap", "Time", "A_c"])
            if write_header:
                writer.writeheader()
            result = {"Order": order_number,
                      "Scenario": scenario_id,
                      "n_s_max": n_s_max,
                      "Solution": round(model.ObjVal / 1000000),
                      "Gap": round(model.MIPGap, 2),
                      "Time": round(model.Runtime),
                      "A_c": round(lower_bound / 1000000)}
            writer.writerow(result)

    else:
        print("No solution found.")
    return model, lower_bound
