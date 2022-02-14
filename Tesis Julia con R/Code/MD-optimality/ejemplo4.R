library(BsMD2)
library(BsMD)
library(JuliaCall)
julia_setup(JULIA_HOME = "C:/Users/Valeria/AppData/Local/Programs/Julia-1.6.3/bin")

setwd("~/ITAM/Tesis/Julia con R/Code/MD-optimality")
julia_source("MDopt.jl")

data(M96e2)
print(M96e2)

set.seed(99)

X <- as.matrix(cbind(blk = rep(-1, 8), M96e2[c(25,2,19,12,13,22,7,32), 1:5]))
y <- M96e2[c(25,2,19,12,13,22,7,32), 6]

pp <- BsProb1(X = X[, 2:6], y = y, p = .25, gamma = .4, 
              max_int = 3, max_fac = 5, top = 32)

p <- pp@p_mod
facs <- pp@fac_mod
Xcand <- as.matrix(cbind(blk = rep(+1, 32), M96e2[, 1:5]))

# Conversiones para los tipos de Julia
X <- as.data.frame(X)
julia_assign("X", X)
julia_assign("y", y)
julia_assign("p_mod", p)
julia_assign("fac_mod", facs)
julia_command("fac_mod = NamedArray(fac_mod)")
julia_eval("fac_mod = Int64.(fac_mod)")
julia_assign("Xcand", Xcand)
julia_command("Xcand = NamedArray(Xcand)")
julia_eval("Xcand = Int64.(Xcand)")

julia_eval("MDopt(X = X, y = y, Xcand = Xcand, nMod = 32, p_mod = p_mod, 
    fac_mod = fac_mod, nFDes = 4, max_int = 3, g = 0.4, Iter = 10, nStart = 25, top = 5)")


