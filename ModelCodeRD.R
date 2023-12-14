##********************************************************************************
## Model code for fitting a rate-distortion model to an empirical confusion matrix
##********************************************************************************
##
## Author: Chris R. Sims (simsc3@rpi.edu)
## Last modified: 2018-02-06
##
## Description: This code accompanies the manuscript 'Efficient Coding Explains the 
## Universal Law of Generalization in Human Perception', and implements fitting
## a model based on rate-distortion theory to an empirical confusion matrix.
## The implements the Blahut algorithm for computing an optimal, but capacity
## limited information channel.
##
## Usage: Calling the function Run() will fit the model to each of the datasets
## listed in the vector all.files. Note, this will take approximately 30 min to
## 1hr per dataset. Run time can be reduced by running multiple instances in
## parallel.
##
## The output of the code is a saved .RDS (R data structure) file containing
## an object with the following variables:
##
## date: The date the model file was created
## model: A string identifying the model
## conf: The empirical confusion matrix to which it was fit
## rho: The estimated cost function
## log.lik: The log likelihood of the model fit
## log.posterior: The log posterior model score
## optim.results: Raw output of the optimization function
## hessian: The computed Hessian matrix
## mu: The mean of the Laplace approximation to the posterior distribution
## sigma: The covariance matrix of the Laplace approximation to the posterior

##********************************************************************************
## Data handling
##********************************************************************************

LoadData <- function(filename) {
    data <- as.matrix(read.csv(filename, header = FALSE))
    colnames(data) <- NULL
    data
}

all.files <- c("Ashby-01",
               "Azadi-Combined",
               "Forrest-Veterinarians",
               "Grey",
               "Hettinger-All",
               "Kornbrot-Biased-Combined",
               "Kornbrot-Unbiased-Combined",
               "Miller-00db",
               "Nosofsky-Color-03",
               "Rouder-13-Monique",
               "Rouder-20-Monique")

##********************************************************************************
## Scripts
##********************************************************************************

Run <- function(safe.mode = TRUE, verbose = FALSE) {

    models <- c("Full")
    attempts <- 3
    
    for(file in all.files) {
        cat("**********", file, "**********\n")
        for(model in models) {
            cat("  ", model, "\n    ")
            input.file <- paste("Data-", file, ".csv", sep = "")
            output.file <- paste("Fit_Bayes_", model, "_", file, ".rds", sep = "")
            
            if(safe.mode) {
                if(file.exists(output.file)) {
                    cat("\n")
                    next
                }
                saveRDS(NULL, output.file)
            }
            conf <- LoadData(input.file)

            opt.val <- -Inf
            for(attempt in 1:attempts) {
                val <- FitModel(conf, model, output.file, verbose = verbose,
                                criterion = opt.val)
                cat(val, " ")
                if(val > opt.val) {
                    opt.val <- val
                }
            }
            cat("\n")
            
            ComputeHessian(file, model)
        }
    }
    return(TRUE)
}

##********************************************************************************
## Model definitions
##********************************************************************************

OptimalChannel <- function(rho,  iters = 100) {
    ## Compute an optimal information channel for a given cost matrix rho.
    ## This function implements the Blahut algorithm (Blahut, 1972).
    
    n <- nrow(rho)
    x <- y <- 1:n

    q <- matrix(data = NA, nrow = n, ncol = n)
    qy <- vector(mode = "numeric", length = n)
    q.new <- matrix(data = NA, nrow = n, ncol = n)
    qy.new <- vector(mode = "numeric", length = n)
    for(i in 1:n) {
        q[i, ] <- rep(1/n, n)
        qy[i] <- 1/n
    }
    
    for(iter in 1:iters) {
        qy.new <- (1/n) * colSums(q)
        for(i in 1:n) {
            q.new[i, ] <- qy * exp(-1 * rho[i, ])
            q.new[i, ] <- q.new[i, ] / sum(q.new[i, ])
        }
        q <- q.new
        qy <- qy.new
    }
    q
}

DataLikelihood <- function(conf, q) {
    ## Compute the log likelihood of an nxn confusion matrix, according to
    ## a conditional probability distribution specified by the nxn matrix q.
    
    n <- nrow(conf)
    log.lik <- 0
    for(i in 1:n) {
        prob <- pmax(1e-10, q[i, ])
        log.lik <- log.lik + dmultinom(conf[i, ], , prob, log = TRUE)
    }
    log.lik
}

MutualInformation <- function(c) {
    ## Computes the mutual information of a confusion matrix
    ## Note: This assumes that stimuli are uniformly distributed
    
    m <- nrow(c);
    n <- ncol(c);
    R <- 0.0;
    px <- vector(mode = "numeric", length = n)
    for(i in 1:n) {
        px[i] <- sum(c[i,])
        c[i, ] <- c[i, ] / max(1, sum(c[i, ]))
    }
    px <- px / sum(px)

    qy <- vector(mode = "numeric", length = n)
    for(j in 1:n) {
        qy[j] <- 0
        for(i in 1:m) {
            qy[j] <- qy[j] + c[i,j] * px[i]
        }
    }
    
    for(i in 1:m) {
        for(j in 1:n) {
            a <- (c[i,j] * px[i])
            if(a > 0) {
                R <- R + a * (log2(c[i,j]) - log2(qy[j]))
            }
        }
    }
    return(R);
}

FindSlope <- function(rho, rate, interval = c(-1e2, 0)) {
    ## For a given cost function and information rate, find the predicted
    ## slope of the generalization gradient
    
    obj.fn <- function(s_est) {
        q <- OptimalChannel(-s_est * rho)
        MutualInformation(q) - rate
    }
    result <- tryCatch({
        uniroot(obj.fn, interval)
    }, error = function(e) {
        NA
    })

    if(is.na(result[1])) {
        return(NA)
    } else {
        result$root
    }
}

##********************************************************************************
## Model fitting
##********************************************************************************

ModelPrior <- function(rho) {
    ## Defines the (log) prior distirbution over the cost function rho
    
    rho.sym <- (1/2)*(rho + t(rho))[upper.tri(rho)]
    rho.asym <- (1/2)*(rho - t(rho))[upper.tri(rho)]
    rho.diag <- diag(rho)
    
    log.prior <-
        sum(dexp(rho.sym, rate = 0.1, log = TRUE)) +
        sum(dexp(rho.diag, rate = 10, log = TRUE)) + 
        sum(dexp(abs(rho.asym), rate = 10, log = TRUE))
    return(log.prior)
}

FitModel <- function(conf, model,
                     output.file = NULL, verbose = FALSE, criterion = -Inf) {
    ## Fit the rate-distortion model to an empirical confusion matrix
    
    n <- nrow(conf)
    rho0 <- matrix(data = rnorm(n*n, mean = 1, sd = 0.1),
                   nrow = n, ncol = n)
    diag(rho0) <- 0
    
    if(model == "Full") {
        par0 <- rho0
    }
    
    obj.fn <- function(par) {
        rho <- ExtractParameters(par, model, n)
        q <- OptimalChannel(rho)
        if(is.na(q[1])) {
            return(-1e10)
        }
        log.prior <- ModelPrior(rho)
        log.lik <- DataLikelihood(conf, q)
        
        log.lik + log.prior
    }
    
    results <- optim(par0, obj.fn,
                     method = "BFGS",
                     control = list(maxit = 10000,
                                    fnscale = -1))
    
    par <- results$par
    rho <- ExtractParameters(par, model, n)
    q <- OptimalChannel(rho)
    log.lik <- DataLikelihood(conf, q)
    log.posterior <- results$value
    
    results <- list(
        date = date(),
        model = model,
        conf = conf,
        rho = rho,
        log.lik = log.lik,
        log.posterior = log.posterior,
        optim.results = results)
    
    if(!is.null(output.file)) {
        if(log.posterior > criterion) {
            saveRDS(results, file = output.file)
        }
    }
    invisible(log.posterior)
}

ExtractParameters <- function(par, model, n) {
    ## Translate a raw parameter vector into an nxn cost matrix
    if(model == "Full") {
        rho <- matrix(data = par, nrow = n, ncol = n)
        rho <- rho - min(rho)
        return(rho)
    }
}

##********************************************************************************
## Misc
##********************************************************************************

ComputeHessian <- function(file, model) {
    ## Computes the Hessian matrix for a given results file
    
    library(numDeriv)
    
    obj.fn <- function(par, conf = NULL, n = NULL) {
        rho <- ExtractParameters(par, model, n)
        q <- OptimalChannel(rho)
        if(is.na(q[1])) {
            return(-1e10)
        }
        log.prior <- ModelPrior(rho)
        log.lik <- DataLikelihood(conf, q)
        
        log.lik + log.prior
    }
    
    input.file <- paste("Fit_Bayes_", model, "_",  file, ".rds", sep = "")
    if(!file.exists(input.file)) {
        cat("File doesn't exist.")
        return(FALSE)
    }
    results <- readRDS(input.file)
    if(is.null(results$rho)) {
        cat("No cost matrix found.")
        return(FALSE)
    }
    if(!is.null(results$hessian)) {
        cat("Hessian already computed.")
        return(FALSE)
    }
    
    n <- nrow(results$conf)
    hess <- hessian(obj.fn, results$optim.results$par,
                    conf = results$conf, n = n)
    
    results$hessian <- hess
    sigma <- solve(-hess)
    results$mu <- as.vector(results$optim.results$par)
    results$sigma <- sigma
    
    saveRDS(results, input.file)
    return(TRUE)
}
