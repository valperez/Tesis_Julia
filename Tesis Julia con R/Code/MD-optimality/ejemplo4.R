# Para guardar los tiempos
tiempos_df <- data.frame(matrix(nrow = 5, ncol = 4))
colnames(tiempos_df) <- c("BsMD2", "BsMD", "JuliaCall", "reticulate")
row.names(tiempos_df) <- seq(1, 5)

runs <- seq(1, 6)


library(BsMD2)

setwd("~/ITAM/Tesis/Julia con R/Code/MD-optimality")
data(M96e2)
print(M96e2)

X <- as.matrix(cbind(blk = rep(-1, 8), M96e2[c(25,2,19,12,13,22,7,32), 1:5]))
y <- M96e2[c(25,2,19,12,13,22,7,32), 6]

pp <- BsProb1(X = X[, 2:6], y = y, p = .25, gamma = .4, 
              max_int = 3, max_fac = 5, top = 32)

p <- pp@p_mod
facs <- pp@fac_mod
Xcand <- as.matrix(cbind(blk = rep(+1, 32), M96e2[, 1:5]))

times <- c()
for (i in runs){
  t <- Sys.time()
  e4_R <- BsMD2::MDopt(X = X, y = y, Xcand = Xcand, 
                       nMod = 32, p_mod = p, fac_mod = facs, 
                       g = 0.4, Iter = 10, nStart = 25, top = 5)
  t_2 <- Sys.time()
  t_final <- difftime(t_2, t, unit = "secs")
  
  times <- c(times, t_final)
}
tiempos_df["BsMD2"] <- times


# # # R paquete original 
library(BsMD)
reactor8.BsProb <- BsProb(X = X, y = y, blk = 1, mFac = 5, mInt = 3, 
                          p = 0.25, g = 0.40, ng = 1, nMod = 32)

nf <- reactor8.BsProb$nftop
s2 <- reactor8.BsProb$sigtop

times <- c()
for (i in runs){
  t_RO <- Sys.time()
  ej4_RO <- BsMD::MD(X = X, y = y, nFac = 5, nBlk = 1, mInt = 3, 
                     g = 0.40, nMod = 32, p = p, s2 = s2, nf = nf, 
                     facs = facs, nFDes = 4, Xcand = Xcand, 
                     mIter = 20, nStart = 25, top = 5)
  t_2 <- Sys.time() 
  t_final <- difftime(t_2, t_RO, unit = "secs")
  times <- c(times, t_final)
}
tiempos_df["BsMD"] <- times


library(JuliaCall)
julia_setup(JULIA_HOME = "C:/Users/Valeria/AppData/Local/Programs/Julia-1.6.3/bin")

julia_source("MDopt.jl")

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

times <- c()
for (i in runs){
  t <- Sys.time()
  julia_eval("MDopt(X = X, y = y, Xcand = Xcand, nMod = 32, p_mod = p_mod, 
    fac_mod = fac_mod, nFDes = 4, max_int = 3, g = 0.4, Iter = 10, nStart = 25, top = 5)")
  t_2 <- Sys.time()
  t_final <- difftime(t_2, t, unit = "secs")
  
  times <- c(times, t_final)
}
tiempos_df["JuliaCall"] <- times

# # # Python con R
library(reticulate)

source_python("MD_Python.py")

X_P <- as.data.frame(X)
Xcand_P <- as.data.frame(Xcand)
fac_mod_P <- as.data.frame(facs)

X_P <- r_to_py(X_P)
y_P <- r_to_py(y) 
Xcand_P <- r_to_py(Xcand_P)
p_mod_P <- r_to_py(p)
fac_mod_P <- r_to_py(fac_mod_P)

nMod_P <- r_to_py(32L)
nFDes_P <- r_to_py(4L)
max_int_P <- r_to_py(3L)
g_P <- r_to_py(0.4)
Iter_P <- r_to_py(10L)
nStart_P <- r_to_py(25L)
top_P <- r_to_py(5L)

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

write.csv(tiempos_df, "tiempos_MD_ej4.csv")


