# Julia con R
# aqui me dice lo de las funciones : https://syl1.gitbook.io/julia-language-a-concise-tutorial/language-core/interfacing-julia-with-other-languages

library(BsMD2)
library(JuliaCall)
julia_setup(JULIA_HOME = "C:/Users/Valeria/AppData/Local/Programs/Julia-1.6.3/bin")

setwd("~/ITAM/Tesis/Julia con R/Code/MD-optimality")
julia_source("MDopt.jl")

set.seed(99)

X <- as.matrix(BM93e3[1:16,c(1,2,4,6,9)]) #matriz de diseÃ±o inicial
y <- as.vector(BM93e3[1:16,10]) #vector de respuesta
p_mod <- c(0.2356,0.2356,0.2356,0.2356,0.0566) #probabilidad posterior de los 5 modelos

fac_mod <- matrix(c(2,1,1,1,1,3,3,2,2,2,4,4,3,4,3,0,0,0,0,4),nrow=5,
                  dimnames=list(1:5,c("f1","f2","f3","f4")))

Xcand <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                  -1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,
                  -1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,
                  -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,
                  -1,1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,-1),
                nrow=16,dimnames=list(1:16,c("blk","f1","f2","f3","f4"))
)

# Conversiones para los tipos de Julia
X <- as.data.frame(X)
julia_assign("X", X)
julia_assign("y", y)
julia_assign("p_mod", p_mod)
julia_assign("fac_mod", fac_mod)
julia_command("fac_mod = NamedArray(fac_mod)")
julia_eval("fac_mod = Int64.(fac_mod)")
julia_assign("Xcand", Xcand)
julia_command("Xcand = NamedArray(Xcand)")
julia_eval("Xcand = Int64.(Xcand)")

#julia_command("X = convert(DataFrame, X)")

julia_eval("MDopt(X = X, y = y, Xcand = Xcand, nMod = 5, 
    p_mod = p_mod, fac_mod = fac_mod, nFDes = 4, max_int = 3, g = 2, Iter = 20, nStart = 10, top = 10)")




