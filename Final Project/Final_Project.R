# Experiment one

library(sandwich)
# Define parameters

set.seed(2811) # For reproducibility

num_simulations <- 1000 # Number of Monte Carlo simulations
sample_size <- c(30, 50, 100, 200, 500) # Number of samples per simulation
# the degree of x used to generate noisy
sample_degree <- c(0, 1 / 5, 1 / 4, 1 / 
                     3, 1 / 2,  1, 2, 3, 4, 5)
# true parameters
beta_0 <- 3
beta_1 <- 5

#' Function to generate data
#'
#' @param n number of observations, x_poly the degree of poly to generate noisy
#' @return the data set contains x, y
generate_data <- function(n, x_poly) {
  x <- runif(n, min=0, max=1)
  noisy <- rnorm(n, mean=0, sd=1)  # Generate n samples from N(0, 1)
  noisy <- noisy * (x ^ x_poly)  # Scale the data by the x_poly
  data_set <- data.frame(x=x, y=beta_0 + beta_1 * x + noisy)
  return(data_set)
}

# MC simulation
# the coverage using normal way
normal_covarage <- matrix(nrow=length(sample_size), ncol=length(sample_degree))
rownames(normal_covarage) <- sample_size
colnames(normal_covarage) <- sample_degree

# the coverage using robust standard error
robust_covarage <- matrix(nrow=length(sample_size), ncol=length(sample_degree))
rownames(robust_covarage) <- sample_size
colnames(robust_covarage) <- sample_degree

for (n in 1: length(sample_size)){
  for (q in 1: length(sample_degree)){
    
    normal_covrage_condition <- c()
    robust_covrage_condition <- c()
    
    for (i in 1: num_simulations){
      df <- generate_data(sample_size[n], sample_degree[q])
      model <- lm(y ~ x, data=df)
      
      # coverage for normal way 
      normal_CI <- confint(model, 'x', level=0.95)
      normal_covrage_condition[i] <- (beta_1 >= normal_CI[1] & beta_1 <= normal_CI[2])
      
      # coverage for Robust standard error
      vhat_HC2 <- vcovHC(model, type="HC2")
      se_HC2 <- sqrt(diag(vhat_HC2))
      t_criticcal <- qt(0.975, sample_size[n] - 2)
      lower_bd <- model$coefficients[2] - t_criticcal * se_HC2[2]
      upper_bd <- model$coefficients[2] + t_criticcal * se_HC2[2]
      robust_CI <- cbind(lower_bd, upper_bd)
      robust_covrage_condition[i] <- (beta_1 >= robust_CI[1] & beta_1 <= robust_CI[2])
    }
    
    normal_covarage[n, q] <- mean(normal_covrage_condition)
    robust_covarage[n, q] <- mean(robust_covrage_condition)
    
  }
} 

covarage_rate <- normal_covarage / robust_covarage
