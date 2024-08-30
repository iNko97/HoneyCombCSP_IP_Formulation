import gurobipy as gp
from gurobipy import GRB
import numpy as np
from order import Order, a_ic_generator

# Initialize model
path = "Input_data.ods"
order_number = 7

model = gp.Model("2D_Cutting_Stock")

# Factory settings
order = Order(path, order_number)

L_min = order.L_min  # Minimum panel length
L_max = order.L_max  # Maximum panel length
n_s_max = order.n_s_max  # Maximum number of stock sizes
n_w_max = order.n_w_max  # Maximum number of widths
n_i_max = order.n_i_max  # Maximum number of items per pattern
one_group = order.one_group  # Only one-groups are allowed
W = order.available_widths  # Set of w available stock widths
I = order.Items  # I: item types with their Width, Length, and Demand


# SETS AND PARAMETERS

# Potential stocks J of dimensions |W| * n_s_max
# This is essentially a matrix of all x_j
# Each row is indexed by available widths.
# Note how for a certain w, idx = W.index(w),then the set J^w' = J[idx, :]
J = np.zeros((len(W), n_s_max))

# Subsets of I compatible with stock size j (C_j)
# Although referring to the J matrix, it only depends on W values.
# given a j at index J[idx, 0], ..., J[idx, n] --> C_j[idx] forall n
C_j = order.C_j_generator()

a_ic = [a_ic_generator(i, c) for i in range(len(I)) for c in range(2 ** len(I) + 1)]
a_ic = np.reshape(a_ic, (len(I), 2**len(I)+1))

# dictionary with triplet (c, idx_j, k) where c I_c, idx_j J.shape[0], k stock size
lmin_cjk = {}
# dictionary with tuple (c, idx_j) : list(range(kmin_cj, ..., kmax_cj))
K_cj = {}

#  dictionary with (idx) : Delta_j, that is the minimum difference between two lmin_cjk
Delta_j = {}

for idx in range(J.shape[0]):
    for c in C_j[idx]:
        (_kmin_cj, _kmax_cj, _lmin_cjk) = order.optimised_stocksize_variables(c, W[idx])

        lmin_cjk.update(_lmin_cjk)
        K_cj[(c, idx)] = list(range(_kmin_cj, _kmax_cj+1))

# PRE-PROCESSING AND MODEL DIMENSIONS REDUCTION

# Proposition 1: Generalised not just for k-1 but for the biggest k' smaller than k
for idx in range(J.shape[0]):
    for c in C_j[idx]:
        for k in K_cj[c, idx]:
            if k == min(K_cj[c, idx]):
                continue
            k_prev = max((x for x in K_cj[c, idx] if x < k))
            if lmin_cjk[c, idx, k_prev] == lmin_cjk[c, idx, k]:
                K_cj[c, idx].remove(k)

# Proposition 2: Substituting I_c with |I_c| cutting pattern with lower lower bounds
for idx in range(J.shape[0]):
    for c in C_j[idx]:
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
                K_cj[c, idx].remove(k)


# \Delta_j, the minimum length difference between two different values of lmin_cjk for a stock j
for idx in range(J.shape[0]):
    best_Delta = 32767
    for c_1 in C_j[idx]:
        for c_2 in C_j[idx]:
            if c_1 & c_2 > 0:
                continue
            K_cj_1 = K_cj[(c_1, idx)]
            K_cj_2 = K_cj[(c_2, idx)]
            for k_1 in K_cj_1:
                for k_2 in K_cj_2:
                    diff = abs(lmin_cjk[(c_1, idx, k_1)] - lmin_cjk[(c_2, idx, k_2)])
                    if diff < best_Delta:
                        best_Delta = diff
                        if best_Delta == 0:
                            break
                else:
                    continue
                break
            else:
                continue
            break
        else:
            continue
        break
    Delta_j[idx] = best_Delta


# DECISION VARIABLES

# \alpha_{cj} 1 if the subset of item types I_c is assigned to stock size j, 0 otherwise for j \in J, c \in C
alpha_cj = model.addVars(
    [(c, idx, n)
     for idx in range(J.shape[0])
     for n in range(J.shape[1])
     for c in C_j[idx]
     ],
    vtype=GRB.BINARY,
    name="alpha_cj"
)

# \beta_j 1 if stock size j is selected, 0 otherwise for j \in J
beta_j = model.addVars(
    [(idx, n)
     for idx in range(J.shape[0])
     for n in range(J.shape[1])
     ],
    vtype=GRB.BINARY,
    name="beta_j"
)

# \gamma_{cjk} 1 if a total of k panels of stock size j are produced to satisfy the demand for item types in subset I_c,
# 0 otherwise, for j \in J, c \in C_j, k \in K_{cj}
gamma_cjk = model.addVars(
    [(c, idx, n, k)
     for idx in range(J.shape[0])
     for n in range(J.shape[1])
     for c in C_j[idx]
     for k in K_cj[(c, idx)]
     ],
    vtype=GRB.BINARY,
    name="gamma_cjk"
)

# \delta_{w} 1 if at least one stock size in J^{w} is selected, 0 otherwise for w' \in W
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
     for idx in range(J.shape[0])
     for n in range(J.shape[1])
     ],
    vtype=GRB.CONTINUOUS,
    lb=0,
    name="x_j"
)

# y_{cjk}
# Auxiliary continuous variable needed to linearise the product x_j*\gamma_{cjk} for j \in J, c \in C, k \in K_{cj}.
# Note that y_{cjk} = x_j if \gamma_{cjk} = 1, otherwise 0
y_cjk = model.addVars(
    [(c, idx, n, k)
     for idx in range(J.shape[0])
     for n in range(J.shape[1])
     for c in C_j[idx]
     for k in K_cj[(c, idx)]],
    vtype=GRB.CONTINUOUS,
    lb=0,
    name="y_cjk"
)

# OBJECTIVE FUNCTION
# 2. Minimise the total area of the material used
model.setObjective(
    gp.quicksum(
        W[idx] * k * y_cjk[c, idx, n, k]
        for idx in range(J.shape[0])
        for n in range(J.shape[1])
        for c in C_j[idx]
        for k in K_cj[(c, idx)]
    ),
    GRB.MINIMIZE
)

# CONSTRAINTS

# 3. Ensure that each item is allocated to a stock size
for i in range(len(I)):
    # noinspection PyTypeChecker
    model.addConstr(
        gp.quicksum(alpha_cj[c, idx, n] * a_ic[i][c]
                    for idx in range(J.shape[0])
                    for n in range(J.shape[1])
                    for c in C_j[idx]) == 1,
        name=f'link_alpha_a_{i}'
    )

# 4. Ensure that a stock size is used iff at least 1 item type is assigned to it
for idx in range(J.shape[0]):
    for n in range(J.shape[1]):
        for c in C_j[idx]:
            model.addConstr(
                alpha_cj[c, idx, n] <= beta_j[idx, n],
                name=f"link_alpha_beta_{c}_{idx}_{n}"
            )

# 5. Ensure the limit on stock sizes is not broken
model.addConstr(
    gp.quicksum(
        beta_j[idx, n]
        for idx in range(J.shape[0])
        for n in range(J.shape[1])
    ) <= n_s_max,
    name="limit_stock_sizes"
)

# 6. Ensure that if a subset I_c is assigned to j, then a k is assigned.
for idx in range(J.shape[0]):
    for n in range(J.shape[1]):
        for c in C_j[idx]:
            model.addConstr(
                gp.quicksum(
                    gamma_cjk[c, idx, n, k]
                    for k in K_cj[(c, idx)]
                    ) == alpha_cj[c, idx, n],
                name=f"link_gamma_alpha_{c}_{idx}_{n}"
            )

# 7. Ensure y_cjk and x_j constraints
# 8. Ensure that the y_cjk respects its length st.
# 9. Ensure that y_cjk doesn't exceed maximum length
# 10. Ensure that x_j doesn't exceed maximum length if k panels are being produced.
for idx in range(J.shape[0]):
    for n in range(J.shape[1]):
        for c in C_j[idx]:
            model.addConstr(
                gp.quicksum(
                    lmin_cjk[(c, idx, k)] * gamma_cjk[c, idx, n, k]
                    for k in K_cj[(c, idx)]
                ) <= x_j[idx, n],
                name=f"length_bound_7_{c}_{idx}_{n}"
            )
            for k in K_cj[(c, idx)]:
                model.addConstr(
                    y_cjk[c, idx, n, k] <= x_j[idx, n],
                    name=f"link_8_y_x_{c}_{idx}_{n}_{k}"
                )
                model.addConstr(
                    y_cjk[c, idx, n, k] <= L_max * gamma_cjk[c, idx, n, k],
                    name=f"link_9_y_x_{c}_{idx}_{n}_{k}"
                )
                model.addConstr(
                    x_j[idx, n] - y_cjk[c, idx, n, k] <= L_max * (1 - gamma_cjk[c, idx, n, k]),
                    name=f"link_10_y_x_{c}_{idx}_{n}_{k}"
                )

# 11. Ensure x_j doesn't exceed maximum length for all selected stock sizes j, otherwise 0.
for idx in range(J.shape[0]):
    for n in range(J.shape[1]):
        model.addConstr(
            L_min * beta_j[idx, n] <= x_j[idx, n],
            name=f"stock_length_min_{idx}_{n}"
        )
        model.addConstr(
            x_j[idx, n] <= L_max * beta_j[idx, n],
            name=f"stock_length_max_{idx}_{n}"
        )

# 12. Ensure a stock size width is used iff used at least one time.
for idx in range(J.shape[0]):
    model.addConstr(
        gp.quicksum(beta_j[idx, n] for n in range(J.shape[1])) >= delta_w[idx],
        name=f"link_delta_beta_{idx}"
    )
    for n in range(J.shape[1]):
        model.addConstr(
            beta_j[idx, n] <= delta_w[idx],
            name=f"link_beta_delta_{idx}_{n}"
        )

# 13. Limit the number of stock widths
model.addConstr(
    gp.quicksum(delta_w[idx] for idx in range(J.shape[0])) <= n_w_max,
    name="limit_stock_widths"
)

# 21. Symmetry breaking constraint for \beta_j
for idx in range(J.shape[0]):
    for n in range(0, J.shape[1]-1):
        model.addConstr(
                    beta_j[idx, n+1] <= beta_j[idx, n],
                    name=f"beta_symmetry_{idx}_{n}"
                )

# 22. Symmetry breaking constraint for \x_j
for idx in range(J.shape[0]):
    for n in range(J.shape[1] - 1):
        model.addConstr(
                    x_j[idx, n] - x_j[idx, n+1] >= Delta_j[idx] * beta_j[idx, n+1],
                    name=f"x_symmetry_{idx}_{n}"
                )


# Optimize the model
model.optimize()

# Extract the solution
if model.status == GRB.OPTIMAL:
    print("Optimal solution found:")
    for idx in range(J.shape[0]):
        for n in range(J.shape[1]):
            if beta_j[idx, n].x > 0.5:
                print(f"Stock size {idx}, {n}: width = {W[idx]}, length = {x_j[idx, n].x}")
                for c in C_j[idx]:
                    if alpha_cj[c, idx, n].x > 0.5:
                        print(f"  Subset {format(c, f'0{len(I)}b')}:")
                        components = [i for i in range(len(I)) if c & (1 << (len(I) - 1 - i))]
                        print(f"  or {components}:")
                        for k in K_cj[(c, idx)]:
                            if gamma_cjk[c, idx, n, k].x > 0.5:
                                print(f"    {k} panels, y_cjk = {y_cjk[c, idx, n, k].x}")
else:
    print("No optimal solution found.")
