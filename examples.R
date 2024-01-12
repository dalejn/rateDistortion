install.packages("RateDistortion")
library(RateDistortion)
install.packages("numDeriv")
library(numDeriv)

# Define a discretized Gaussian information source
x <- seq(from = -10, to = 10, length.out = 100)
Px <- dnorm(x, mean = 0, sd = 3)
Px <- Px / sum(Px) # Ensure that probability sums to 1
y <- x # The destination alphabet is the same as the source

# Blahut Algorithm --------------------------------------------------------

# Define a quadratic cost function
cost.function <- function(x, y) {
  (y - x)^2
}

# A different cost function with additional named input arguments
cost.function.2 <- function(x, y, alpha = 1.0) {
  abs(y - x) ^ alpha
}

# Slope of the rate-distortion curve
s <- -1

# Compute the rate-distortion value at the given point s
channel <- BlahutAlgorithm(x, Px, y, cost.function, s)

# Not run:
# channel.2 <- BlahutAlgorithm(x, Px, y, cost.function.2, s, alpha = 1.5)

print(channel)


# Channel distortion ------------------------------------------------------

# Compute the channel distortion according to a different cost function
abs.cost.function <- function(x, y) {
  abs(y - x)
}

ChannelDistortion(channel, abs.cost.function)


# Conditional distribution ------------------------------------------------

# Compute the rate-distortion value at the given point s
channel <- BlahutAlgorithm(x, Px, y, cost.function, s)

# Compute & plot the conditional probability distribution for a particular channel input.
cpd <- ConditionalDistribution(channel, index = 50)
plot(cpd$y, cpd$p)


# Channel difference (or error) distribution ------------------------------

# Compute & plot the channel difference (or error) distribution
diff.distribution <- DifferenceDistribution(channel)
plot(diff.distribution$diff, diff.distribution$p)

# Find optimal channel ----------------------------------------------------

R <- 1.5
Q <- FindOptimalChannel(x, Px, y, cost.function, R, verbose = TRUE)
print(Q)

# Find rate ---------------------------------------------------------------

D <- 1.0
Q <- FindRate(x, Px, y, cost.function, D, verbose = TRUE)
print(Q)


# Maximum cost ------------------------------------------------------------

MaximumCost(x, y, cost.function)


# Mutual information ------------------------------------------------------

# Compute the information rate for this channel assuming a different (uniform) source distribution
uniform.dist <- rep(1 / 100, 100)
MutualInformation(channel, uniform.dist)


# Draw random samples of the conditional probability distribution --------

# Draw random samples from the output of this channel, for a given input
samples <- Sample(channel, 1000, 50, show.progress = TRUE)

# Plot a histogram of the output
hist(samples$y)
