# Para guardar los tiempos
tiempos_df <- data.frame(matrix(nrow = 5, ncol = 4))
colnames(tiempos_df) <- c("BsMD2", "BsMD", "JuliaCall", "reticulate")
row.names(tiempos_df) <- seq(1, 5)

runs <- seq(1, 5)

# # # R tesis Paty
library(BsMD2)
setwd("~/ITAM/Tesis/Julia con R/Code/MD-optimality")

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

times <- c()
for (i in runs){
  t <- Sys.time()
  e3_R <- BsMD2::MDopt(X = X, y = y, Xcand = Xcand,
                       nMod = 5, p_mod = p_mod, fac_mod = fac_mod, 
                       nStart = 25)
  t_2 <- Sys.time()
  t_final <- difftime(t_2, t, unit = "secs")
  
  times <- c(times, t_final)
}
tiempos_df["BsMD2"] <- times


# # # R paquete original
library(BsMD)

s2 <- c(0.5815, 0.5815, 0.5815, 0.5815, 0.4412)

times <- c()
for (i in runs){
  t_RO <- Sys.time()
  e3_RO <- BsMD::MD(X = X, y = y, nFac = 4, nBlk = 1, mInt = 3, g = 2, 
                    nMod = 5, p = p_mod, s2 = s2, nf = c(3, 3, 3, 3, 4), 
                    facs = fac_mod, nFDes = 4, Xcand = Xcand, mIter = 20, 
                    nStart = 25, top = 10)
  t_2 <- Sys.time() 
  t_final <- difftime(t_2, t_RO, unit = "secs")
  times <- c(times, t_final)
}
tiempos_df["BsMD"] <- times


# # #  Julia con R
library(JuliaCall)
julia_setup(JULIA_HOME = "C:/Users/Valeria/AppData/Local/Programs/Julia-1.6.3/bin")

julia_source("MDopt.jl")
# Conversiones para los tipos de Julia
X_J <- as.data.frame(X)
julia_assign("X_J", X_J)
julia_assign("y_J", y)
julia_assign("p_mod_J", p_mod)
julia_assign("fac_mod_J", fac_mod)
julia_command("fac_mod_J = NamedArray(fac_mod_J)")
julia_eval("fac_mod_J = Int64.(fac_mod_J)")
julia_assign("Xcand_J", Xcand)
julia_command("Xcand_J = NamedArray(Xcand_J)")
julia_eval("Xcand_J = Int64.(Xcand_J)")


times <- c()
for (i in runs){
  t <- Sys.time()
  julia_eval("MDopt(X = X_J, y = y_J, Xcand = Xcand_J, nMod = 5, 
    p_mod = p_mod_J, fac_mod = fac_mod_J, nFDes = 4, max_int = 3, g = 2, Iter = 20, nStart = 10, top = 10)")
  t_2 <- Sys.time()
  t_final <- difftime(t_2, t, unit = "secs")
  
  times <- c(times, t_final)
}
tiempos_df["JuliaCall"] <- times



# Python con R
library(reticulate)

source_python("MD_Python.py")

X_P <- as.data.frame(X)
Xcand_P <- as.data.frame(Xcand)
fac_mod_P <- as.data.frame(fac_mod)

X_P <- r_to_py(X_P)
y_P <- r_to_py(y) 
Xcand_P <- r_to_py(Xcand_P)
p_mod_P <- r_to_py(p_mod)
fac_mod_P <- r_to_py(fac_mod_P)

nMod_P <- r_to_py(5L)
nFDes_P <- r_to_py(4L)
max_int_P <- r_to_py(3L)
g_P <- r_to_py(2L)
Iter_P <- r_to_py(20L)
nStart_P <- r_to_py(25L)
top_P <- r_to_py(10L)

times <- c()
for (i in runs){
  t <- Sys.time()
  MD_Python(X = X_P, y = y_P, Xcand = Xcand_P, nMod = nMod_P, 
            p_mod = p_mod_P, fac_mod = fac_mod_P, 
            nFDes = nFDes_P, max_int = max_int_P, 
            g = g_P, Iter = Iter_P, nStart = nStart_P, top = top_P)
  t_2 <- Sys.time()
  t_final <- difftime(t_2, t, unit = "secs")
  
  times <- c(times, t_final)
}
tiempos_df["reticulate"] <- times

write.csv(tiempos_df, "tiempos_MD_ej3.csv")

