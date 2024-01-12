# Function to calculate Shepard's generalization measure from a frequency matrix
shepards_generalization_matrix <- function(frequency_matrix) {
  # Assuming frequency_matrix is a square matrix
  n <- nrow(frequency_matrix)
  probability_matrix <- frequency_matrix / rowSums(frequency_matrix)
  # Initialize a matrix for the generalization measures
  G_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Loop through all pairs of stimuli
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate the generalization measure for each pair
      G_matrix[i, j] <- sqrt((probability_matrix[i, j] * probability_matrix[j, i]) / 
                               (probability_matrix[i, i] * probability_matrix[j, j]))
    }
  }
  
  return(G_matrix)
}

shepards_generalization_matrix(df$conf)
plot(df$rho, shepards_generalization_matrix(df$conf))

#RMSE
actual <- df$conf/rowSums(df$conf)
predicted <- OptimalChannel(df$rho)
sqrt(mean((actual - predicted)^2))
plot(predicted, actual)

empirical <- shepards_generalization_matrix(df$conf)
s <- FindSlope(df$rho, 0.5)
theoretical <- exp(s*df$rho)
sqrt(mean((empirical - theoretical)^2))
plot(theoretical, empirical)

# RD curve
plot(df$rho, theoretical)
