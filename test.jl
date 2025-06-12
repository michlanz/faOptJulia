using Pkg
Pkg.add("JuMP");
Pkg.add("HiGHS");
Pkg.add("DataFrames");
Pkg.add("CSV");
Pkg.add("Plots");

using JuMP, HiGHS, DataFrames, CSV, Plots, Printf;

# import
L = 1 #peso
na = 14 #numero anime
np = 14 #numero famiglie prodotto
nm = 2 #numero materiali
m = 5 #numero macchine

V = 10_000 #ciclo vita degli stampi
I = 600_000 #costo di adattamento di una macchina

C = [3000, 60000] #Cj
Q = Matrix(CSV.read("param_Qij.csv", DataFrame, header = false)); #Qij
U = Matrix(CSV.read("param_Uik.csv", DataFrame, header = false)); #Uik
p = Matrix(CSV.read("param_Pk.csv", DataFrame, header = false))[:,1];  #Pk

P = U * p; #così è Pi
Ptot = sum(P);
M = Ptot / m; #M produzione standard di una macchina

Yub = maximum(P)/V;
Kub = m;

#
# Stesura modello
# Definizione del modello
sandModel = Model(HiGHS.Optimizer)
set_silent(sandModel)

set_optimizer_attribute(sandModel, "mip_rel_gap", 1e-4) # Set relative MIP gap tolerance to a lower value
set_optimizer_attribute(sandModel, "presolve", "on") # Enable presolve
set_time_limit_sec(sandModel, 10.0) # Set time limit for HiGHS


# Definizione delle variabili
@variable(sandModel, X[1:na, 1:2] >= 0)  # Variabile continua
@variable(sandModel, Z[1:na, 1:2], Bin) #attivazione produzione
@variable(sandModel, 0<= Y[1:na, 1:2] <= Yub + 1, Int)  # Variabile intera positiva
@variable(sandModel, K >=0 , Int)  # Variabile intera positiva

# Vincoli
@constraint(sandModel, prod_a[i in 1:na], sum(X[i, j] for j in 1:2) >= P[i])
@constraint(sandModel, prod_b[i in 1:na, j in 1:2], X[i, j] <= P[i] * Z[i, j])
@constraint(sandModel, prodActivation[i in 1:na], sum(Z[i, j] for j in 1:2) .== 1)
@constraint(sandModel, mold[i in 1:na, j in 1:2], X[i, j] <= V * Y[i, j])
@constraint(sandModel, machine, sum(X[i, 2] for i in 1:na) <= K * M )
@constraint(sandModel, pollution, sum(Q[i, j] * X[i, j] for i in 1:na, j in 1:2) <= L * sum(P[i] * Q[i, 1] for i in 1:na) + (1 - L) * sum(P[i] * Q[i, 2] for i in 1:na))

# Funzione obiettivo
@objective(sandModel, Min, sum(C[j] * Y[i, j] for i in 1:na, j in 1:2) + I*K)

vett_y_sm = [] #qui le emissioni
vett_x_sm = [] #qui i soldi spesi
vett_k_sm = [] #qui metto le macchine

for L in 1:-0.05:0.00
  set_normalized_rhs(pollution, L * sum(P[i] * Q[i, 1] for i in 1:na) + (1 - L) * sum(P[i] * Q[i, 2] for i in 1:na))
  optimize!(sandModel)

  push!(vett_x_sm, objective_value(sandModel))
  push!(vett_y_sm, sum(Q[i, j] * value(X[i, j]) for i in 1:na, j in 1:2))
  push!(vett_k_sm, value(K))
  println("Istanza ", L, "      Gap: ", @sprintf("%.0e%%", 100 * relative_gap(sandModel)),"     K = ", value(K))
  #println("MIP gap: ", relative_gap(sandModel))
  #println("K: ", value(K))
  #println(objective_value(sandModel))
  # per bloccare il ciclo quando un modello è infeasible: "termination_status(toyModel) == MOI.INFEASIBLE"
end
println()
println()

#stampare
x_ticks = range(0, stop=100, step=10);#
#
modified_cost = (vett_x_sm ./ minimum(vett_x_sm));
#modified_cost = log10.(vett_x_sm)
plot(modified_cost, 100*vett_y_sm/maximum(vett_y_sm),
    xlabel="Soldi spesi",
    ylabel="Emissioni",
    title="Frontiera di Pareto norm. percentuale",
    label="Pareto",
    legend= (0.85, 0.7),
#    xticks = (x_ticks, x_ticks),
#    yticks = (x_ticks, x_ticks)
)

x_first = modified_cost[1]
y_first = 100*vett_y_sm[1]/maximum(vett_y_sm)
x_last = modified_cost[end]
y_last = 100*vett_y_sm[end]/maximum(vett_y_sm)

# Plot the bisector line
plot!([x_first, x_last], [y_first, y_last],
    label="Bisettrice",
    linestyle=:dash,
    color=:black
)

# Create a secondary y-axis on the right
plot2 = twinx()

# Plot the values of K on the secondary y-axis with magenta color
plot!(plot2, modified_cost, vett_k_sm, # Using vett_k_sm
      label="K",
      ylabel="K (Machines)",
      legend= (0.94, 0.85),
      linestyle=:dashdot,
      color=:magenta) # Set color to magenta