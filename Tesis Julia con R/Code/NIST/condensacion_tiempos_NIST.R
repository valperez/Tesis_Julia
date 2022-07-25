# # # Script para condensar los resultados de los tiempos 
# # # para el ejercicio de NIST
library(polynom)
library(knitr)
library(pracma)
library(dplyr)

setwd("~/GitHub/Tesis_Julia/Tesis Julia con R/Code/NIST")

# Lectura de tiempos de R
r_data <- read.csv("tiempos_NISt_R_Condensado.csv")
r_data <- r_data[, c("Grado", "Promedio")]
colnames(r_data) <- c("Grado", "R")

# Lectura de tiempos de Julia
setwd("~/GitHub/Tesis_Julia/Tesis Julia con R/Code/NIST/Resultados_Julia")
temp <- list.files(pattern = "tiempos_grado")
myfiles <- lapply(temp, read.csv)

julia_data <- data.frame(matrix(NA, nrow = 10, ncol = 4))

for (i in 1:length(myfiles)){
  aux <- myfiles[[i]]
  julia_data[i, ] <- aux[which(aux$X == "mean"), 2:5]
}

nombres_columnas <- colnames(myfiles[[1]])
nombres_columnas <- nombres_columnas[2:length(nombres_columnas)]
colnames(julia_data) <- nombres_columnas 

# Lectura de tiempos Julia
setwd("~/GitHub/Tesis_Julia/Tesis Julia con R/Code/NIST/Python")

python_data <- read.csv("tiempos_NISt_Python.csv")
python_data <- python_data[, c("Grado", "Promedio")]
colnames(python_data) <- c("Grado", "Python")


# Concatenación de todo
tiempos_df <- cbind(r_data[, "Grado"], 
                    julia_data, 
                    r_data[, "R"],
                    python_data[, "Python"])
colnames(tiempos_df) <- c("Grado", 
                          colnames(julia_data),
                          "R", "Python")
                    
