# load libraries
library(dagitty)

# make dags
lag_dag <- dagitty('dag {
  "Ecosystem Function (0)" [pos="-0.468,-0.003"]
  "Ecosystem Function (t)" [pos="0.641,-0.199"]
  "Functional Stability" [outcome,pos="0.084,0.393"]
  "Mean Soil Nutrients" [pos="-0.058,-1.608"]
  "Response Diversity (Temperature)" [exposure,pos="-1.607,0.440"]
  "Species Function (0, 0)" [pos="-0.829,-0.654"]
  "Species Function (0, t)" [pos="0.366,-0.671"]
  "Species Function (n, 0)" [pos="-0.272,-0.667"]
  "Species Function (n, t)" [pos="0.953,-0.654"]
  "Species Performance (0, 0)" [pos="-0.827,-1.127"]
  "Species Performance (0, t)" [pos="0.364,-1.122"]
  "Species Performance (n, 0)" [pos="-0.259,-1.131"]
  "Species Performance (n, t)" [pos="0.959,-1.110"]
  "Species Richness" [adjusted,pos="-1.416,-0.343"]
  "Total Ecosystem Function (AUC)" [pos="0.929,0.406"]
  Asynchrony [pos="-0.723,0.372"]
  "Ecosystem Function (0)" -> "Functional Stability"
  "Ecosystem Function (0)" -> "Total Ecosystem Function (AUC)"
  "Ecosystem Function (t)" -> "Functional Stability"
  "Ecosystem Function (t)" -> "Total Ecosystem Function (AUC)"
  "Mean Soil Nutrients" -> "Species Performance (0, 0)"
  "Mean Soil Nutrients" -> "Species Performance (0, t)"
  "Mean Soil Nutrients" -> "Species Performance (n, 0)"
  "Mean Soil Nutrients" -> "Species Performance (n, t)"
  "Response Diversity (Temperature)" -> Asynchrony
  "Species Function (0, 0)" -> "Ecosystem Function (0)"
  "Species Function (0, t)" -> "Ecosystem Function (t)"
  "Species Function (n, 0)" -> "Ecosystem Function (0)"
  "Species Function (n, t)" -> "Ecosystem Function (t)"
  "Species Performance (0, 0)" -> "Species Function (0, 0)"
  "Species Performance (0, t)" -> "Species Function (0, t)"
  "Species Performance (n, 0)" -> "Species Function (n, 0)"
  "Species Performance (n, t)" -> "Species Function (n, t)"
  "Species Richness" -> "Ecosystem Function (0)"
  "Species Richness" -> "Ecosystem Function (t)"
  "Species Richness" -> "Response Diversity (Temperature)"
  "Species Richness" -> "Species Performance (0, 0)"
  "Species Richness" -> "Species Performance (0, t)"
  "Species Richness" -> "Species Performance (n, 0)"
  "Species Richness" -> "Species Performance (n, t)"
  "Species Richness" -> Asynchrony
  "Total Ecosystem Function (AUC)" -> "Functional Stability"
  Asynchrony -> "Functional Stability"
}')

plot(lag_dag)