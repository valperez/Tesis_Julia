# Vamos a aprender a calcular la descomposición de valores singulares
# de una matriz

# Primero vamos a bajar los paquetes
using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("CSV")
Pkg.add("DataFrames")
using LinearAlgebra, CSV, DataFrames

filip = CSV.read("filip_data.csv", DataFrame)
println("Supongamos que tenemos la siguiente matriz: ")
display(filip)

df = Matrix(filip)

println("Ahora obtengamos la descomposicion de valores singulares de A ")
F = svd(df)
display(F)

println("Vamos a descontruirlo en variables separadas")
u, s, v = F

println("Entonces la matriz U es")
display(u)

println("Los valores singulares son")
display(s)

println("Finalmente, la matriz V es")
display(v)
