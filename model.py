import gurobipy as gp
from gurobipy import GRB
import numpy as np
from preprocessor import C_j_generator, a_ic_generator, optimised_stocksize_variables

# Initialize model
model = gp.Model("2D_Cutting_Stock")

# Factory settings

L_min = 700  # Minimum panel length
L_max = 3100  # Maximum panel length
n_s_max = 2  # Maximum number of stock sizes
n_w_max = 2  # Maximum number of widths
W = [1200, 1400, 1550, 1600]  # Set of available stock widths

# I: item types with their Width, Length, and Demand
I = [
    [230, 250, 600],
    [230, 2140, 600],
    [290, 2140, 600]
]

# Sets and parameters

# Potential stocks J of dimensions |W| * n_s_max
# This is essentially a matrix of all x_j
# Each row is indexed by available widths.
# Note how for a certain w, idx = W.index(w),then the set J^w' = J[idx, :]
J = np.zeros((len(W), n_s_max))

# Subsets of I compatible with stock size j (C_j)
# Although referring to the J matrix, it only depends on W values.
# given a j at index J[idx, 0], ..., J[idx, n] --> C_j[idx] forall n
C_j = C_j_generator()

a_ic = [a_ic_generator(i, c) for i in range(len(I)) for c in range(2**len(I))]
a_ic = np.reshape(a_ic, (len(I), 2**len(I)))

lmin_cjk = {}
kmin_cj = {}
kmax_cj = {}

for idx in range(J.shape[0]):
    for c in C_j[idx]:
        print(W[idx])
        (_kmin_cj, _kmax_cj, _lmin_cjk) = optimised_stocksize_variables(c, W[idx])

        lmin_cjk.update(_lmin_cjk)
        kmin_cj[(c, idx)] = _kmin_cj
        kmax_cj[(c, idx)] = _kmax_cj


# Decision variables
# \alpha_{cj} 1 if the subset of item types I_c is assigned to stock size j, 0 otherwise for j \in J, c \in C
alpha_cj = model.addVars(J, C_j, vtype=GRB.BINARY, name="alpha_cj")

# \beta_j 1 if stock size j is selected, 0 otherwise for j \in J
beta_j = model.addVars(J, vtype=GRB.BINARY, name="beta_j")

# \gamma_{cjk} 1 if a total of k panels of stock size j are produced to satisfy the demand for item types in subset I_c,
# 0 otherwise, for j \in J, c \in C_j, k \in K_{cj}
gamma_cjk = model.addVars(
    [(c, J[idx, n], k)
     for idx in range(J.shape[0])
     for n in range(J.shape[1])
     for c in C_j[idx]
     for k in range(kmin_cj[(c, J[idx, n])], kmax_cj[(c, J[idx, n])] + 1)],
    vtype=GRB.BINARY,
    name="gamma_cjk"
)

# \delta_{w} 1 if at least one stock size in J^{w} is selected, 0 otherwise for w' \in W
delta_w = model.addVars(W, vtype=GRB.BINARY, name="delta_w")

# x_j length of stock size j, for j \in J
x_j = model.addVars(J, vtype=GRB.CONTINUOUS, name="x_j")

# y_{cjk}
# Auxiliary continuous variable needed to linearise the product x_j*\gamma_{cjk} for j \in J, c \in C, k \in K_{cj}.
# Note that y_{cjk} = x_j if \gamma_{cjk} = 1, otherwise 0
y_cjk = model.addVars(
    [(c, J[idx, n], k)
     for idx in range(J.shape[0])
     for n in range(J.shape[1])
     for c in C_j[idx]
     for k in range(kmin_cj[(c, J[idx, n])], kmax_cj[(c, J[idx, n])] + 1)],
    vtype=GRB.CONTINUOUS,
    name="y_cjk"
)

# Objective function
# 2. Minimise the total area of the material used
model.setObjective(
    gp.quicksum(
        W[idx] * k * y_cjk[c, J[idx, n], k]
        for idx in range(J.shape[0])
        for n in range(J.shape[1])
        for c in C_j[idx]
        for k in range(kmin_cj[(c, J[idx, n])], kmax_cj[(c, J[idx, n])] + 1)
    ),
    GRB.MINIMIZE
)

# Constraints
# 3. Ensure that each item is allocated to a stock size
for i in I:
    model.addConstr(gp.quicksum(alpha_cj[c, j] * a_ic[(i, c)] for j in J for c in C_j) == 1, name=f"link_alpha_a_{i}")

# 4. Ensure that a stock size is used iff at least 1 item type is assigned to it
for j in J:
    for c in C_j[j]:
        model.addConstr(alpha_cj[c, j] <= beta_j[j], name=f"link_alpha_beta_{c}_{j}")

# 5. Ensure the limit on stock sizes is not broken
model.addConstr(gp.quicksum(beta_j[j] for j in J) <= n_s_max, name="limit_stock_sizes")

# 6. Ensure that if a subset I_c is assigned to j, then a k is assigned.
for j in J:
    for c in C_j[j]:
        model.addConstr(gp.quicksum(gamma_cjk[c, j, k] for k in range(kmin_cj[(c, j)], kmax_cj[(c, j)] + 1)) == alpha_cj[c, j],
                        name=f"link_gamma_alpha_{c}_{j}")

# 7. Ensure y_cjk and x_j constraints
# 8. Ensure that the y_cjk respects its length st.
# 9. Ensure that y_cjk doesn't exceed maximum length
# 10. Ensure that x_j doesn't exceed maximum length if k panels are being produced.
for j in J:
    for c in C_j[j]:
        model.addConstr(gp.quicksum(lmin_cjk[(c, j, k)] * gamma_cjk[c, j, k] for k in range(kmin_cj[(c, j)], kmax_cj[(c, j)] + 1)) <= x_j[j], name=f"length_bound_{c}_{j}")
        for k in range(kmin_cj[(c, j)], kmax_cj[(c, j)] + 1):
            model.addConstr(y_cjk[c, j, k] <= x_j[j], name=f"link_y_x_1_{c}_{j}_{k}")
            model.addConstr(y_cjk[c, j, k] <= L_max * gamma_cjk[c, j, k], name=f"link_y_x_2_{c}_{j}_{k}")
            model.addConstr(x_j[j] - y_cjk[c, j, k] <= L_max * (1 - gamma_cjk[c, j, k]), name=f"link_y_x_3_{c}_{j}_{k}")

# 11. Ensure x_j doesn't exceed maximum length for all selected stock sizes j, otherwise 0.
for j in J:
    model.addConstr(L_min * beta_j[j] <= x_j[j], name=f"stock_length_min_{j}")
    model.addConstr(x_j[j] <= L_max * beta_j[j], name=f"stock_length_max_{j}")

# 12. Ensure a stock size width is used iff used at least one time.
for w in W:
    model.addConstr(gp.quicksum(beta_j[j] for j in J[w]) >= delta_w[w], name=f"link_delta_beta_{w}")
    for j in J[w]:
        model.addConstr(beta_j[j] <= delta_w[w], name=f"link_beta_delta_{w}_{j}")

# 13. Limit the number of stock widths
model.addConstr(gp.quicksum(delta_w[w] for w in W) <= n_w_max, name="limit_stock_widths")

# Optimize the model
model.optimize()

# Extract the solution
if model.status == GRB.OPTIMAL:
    print("Optimal solution found:")
    for j in J:
        if beta_j[j].x > 0.5:
            print(f"Stock size {j}: length = {x_j[j].x}")
            for c in C_j:
                if alpha_cj[c, j].x > 0.5:
                    print(f"  Subset {c}:")
                    for k in range(kmin_cj[(c, j)], kmax_cj[(c, j)] + 1):
                        if gamma_cjk[c, j, k].x > 0.5:
                            print(f"    {k} panels, y_cjk = {y_cjk[c, j, k].x}")
else:
    print("No optimal solution found.")
