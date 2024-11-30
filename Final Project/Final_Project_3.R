# Experiment two

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
beta_2 <- 4
beta_3 <- 6

#' Function to generate data
#'
#' @param n number of observations, x_poly the degree of poly to generate noisy
#' @return the data set contains x_1, x_2, y
generate_data <- function(n, x_poly) {
  x_1 <- runif(n, min=0, max=1)
  x_2 <- runif(n, min=0, max=1)
  x_3 <- runif(n, min=0, max=1)
  noisy <- rnorm(n, mean=0, sd=1)  # Generate n samples from N(0, 1)
  noisy <- noisy * (x_1 ^ x_poly)  # Scale the data by the x_poly
  data_set <- data.frame(x_1=x_1, x_2=x_2, x_3=x_3, 
                         y=beta_0 + beta_1 * x_1 + beta_2 * x_2 + beta_3 * x_3 + noisy)
  return(data_set)
}

# MC simulation
# the coverage using normal way
normal_covarage_2 <- matrix(nrow=length(sample_size), ncol=length(sample_degree))
rownames(normal_covarage_2) <- sample_size
colnames(normal_covarage_2) <- sample_degree

# the coverage using robust standard error
robust_covarage_2 <- matrix(nrow=length(sample_size), ncol=length(sample_degree))
rownames(robust_covarage_2) <- sample_size
colnames(robust_covarage_2) <- sample_degree

# the coverage using normal way
normal_covarage_3 <- matrix(nrow=length(sample_size), ncol=length(sample_degree))
rownames(normal_covarage_3) <- sample_size
colnames(normal_covarage_3) <- sample_degree

# the coverage using robust standard error
robust_covarage_3 <- matrix(nrow=length(sample_size), ncol=length(sample_degree))
rownames(robust_covarage_3) <- sample_size
colnames(robust_covarage_3) <- sample_degree

for (n in 1: length(sample_size)){
  for (q in 1: length(sample_degree)){
    
    normal_covrage_condition_2 <- c()
    robust_covrage_condition_2 <- c()
    
    normal_covrage_condition_3 <- c()
    robust_covrage_condition_3 <- c()
    
    for (i in 1: num_simulations){
      df <- generate_data(sample_size[n], sample_degree[q])
      model <- lm(y ~ x_1 + x_2 + x_3, data=df)
      
      # coverage for normal way 
      normal_CI_2 <- confint(model, 'x_2', level=0.95)
      normal_covrage_condition_2[i] <- (beta_2 >= normal_CI_2[1] & beta_2 <= normal_CI_2[2])
      
      normal_CI_3 <- confint(model, 'x_3', level=0.95)
      normal_covrage_condition_3[i] <- (beta_3 >= normal_CI_3[1] & beta_3 <= normal_CI_3[2])
      
      # coverage for Robust standard error
      vhat_HC2 <- vcovHC(model, type="HC2")
      se_HC2 <- sqrt(diag(vhat_HC2))
      t_criticcal <- qt(0.975, sample_size[n] - 3)
      
      lower_bd_2 <- model$coefficients[3] - t_criticcal * se_HC2[3]
      upper_bd_2 <- model$coefficients[3] + t_criticcal * se_HC2[3]
      
      lower_bd_3 <- model$coefficients[4] - t_criticcal * se_HC2[4]
      upper_bd_3 <- model$coefficients[4] + t_criticcal * se_HC2[4]
      
      robust_CI_2 <- cbind(lower_bd_2, upper_bd_2)
      robust_covrage_condition_2[i] <- (beta_2 >= robust_CI_2[1] & beta_2 <= robust_CI_2[2])
      
      robust_CI_3 <- cbind(lower_bd_3, upper_bd_3)
      robust_covrage_condition_3[i] <- (beta_3 >= robust_CI_3[1] & beta_3 <= robust_CI_3[2])
    }
    
    normal_covarage_2[n, q] <- mean(normal_covrage_condition_2)
    robust_covarage_2[n, q] <- mean(robust_covrage_condition_2)
    
    normal_covarage_3[n, q] <- mean(normal_covrage_condition_3)
    robust_covarage_3[n, q] <- mean(robust_covrage_condition_3)
    
  }
} 

covarage_rate_2 <- normal_covarage_2 / robust_covarage_2
covarage_rate_3 <- normal_covarage_3 / robust_covarage_3
