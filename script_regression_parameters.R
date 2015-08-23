#clean workspace
rm(list = ls()) 
library(ggplot2) #dataviz
library(gridExtra) #dataviz

#set working directory to the current source file directory
dir.wd <- "/Users/Marco/Dropbox/R workspace/github/regression_parameters/Regression-parameter-estimates/"

#load dataset
df_hrv_age <- read.csv(paste(dir.wd, "df_hrv_age.csv", sep = ""), header = TRUE)

#plot data (scatterplot with distributions)
#placeholder plot - prints nothing at all
margin_plots <- 0.2
empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + 
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank())

# scatterplot of x and y variables
scatter <- ggplot(df_hrv_age, aes(age, mean_hrv)) + geom_point() + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Heart Rate Variability (rMSSD)") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1))

# marginal density of x - plot on top
plot_top <- ggplot(df_hrv_age, aes(age)) + geom_histogram() + 
  scale_x_continuous("Age") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = "none")  

# marginal density of y - plot on the right
plot_right <- ggplot(df_hrv_age, aes(mean_hrv)) + geom_histogram() + 
  coord_flip()  + 
  scale_x_continuous("Heart Rate Variability (rMSSD)") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = "none") 

# arrange the plots together, with appropriate height and width for each row and column
grid.arrange(plot_top, empty, scatter, plot_right, ncol = 2, nrow = 2, 
             widths = c(4, 1), heights = c(1, 4))

#correlation between variables
rcorr(df_hrv_age$age, df_hrv_age$mean_hrv)

#1 - Least squares
#Manual implementation

#add intercept
X <- cbind(rep(1, length(df_hrv_age$age)), df_hrv_age$age)
y <- df_hrv_age$mean_hrv
  
#use correct operator for matrix multiplcation
#also solve can be used because we have a square matrix after multiplying X for its transpose
B_hat <- (solve(t(X)%*%X) %*% t(X)) %*% y
print(B_hat)

#lm package 
lm_coefficients <- summary(lm(mean_hrv ~ age, data = df_hrv_age))$coefficients
print(lm_coefficients)

#plot scatterplot and regression line
p1 <- ggplot(df_hrv_age, aes(age, mean_hrv)) + geom_point() + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Heart Rate Variability (rMSSD)") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  geom_abline(intercept = B_hat[1], slope = B_hat[2], col = "darkred") +
  ggtitle("Manual implementation")
p2 <- ggplot(df_hrv_age, aes(age, mean_hrv)) + geom_point() + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Heart Rate Variability (rMSSD)") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  geom_abline(intercept = lm_coefficients[1], slope = lm_coefficients[2], col = "darkred") +
  ggtitle("lm package")
p1
p2


#2 - Gradient descent
X <- cbind(rep(1, length(df_hrv_age$age)), df_hrv_age$age)
y <- df_hrv_age$mean_hrv

#define cost function
cost <- function(theta_temp)
{
  return (sum(((X%*%theta_temp)-y)^2)/(2*m))
}

#initialize parameters
alpha <- 0.0001
iterations <- 1000000
theta<- c(0, 0)

#samples
m <- nrow(X)

#gradient descent
history <- matrix(NA, iterations, 3) #theta1, theta2, cost
cost_history <- c()
for(i in 1:iterations)
{
  theta[1] <- theta[1] - alpha * (1/m) * sum(((X %*% theta) - y))
  theta[2] <- theta[2] - alpha * (1/m) * sum(((X %*% theta) - y) * X[,2])
  history[i, ] <- c(theta[1], theta[2], cost(theta))
}

p1 <- ggplot(df_hrv_age, aes(age, mean_hrv)) + geom_point() + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Heart Rate Variability (rMSSD)") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  geom_abline(intercept = history[1, 1], slope = history[1, 2], col = "darkred") +
  ggtitle("Gradient descent, first iteration")

p2 <- ggplot(df_hrv_age, aes(age, mean_hrv)) + geom_point() + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Heart Rate Variability (rMSSD)") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  geom_abline(intercept = history[iterations, 1], slope = history[iterations, 2], col = "darkred") +
  ggtitle("Gradient descent, last iteration")
#parameters estimates
p1
p2

#parameters space
df_history <- data.frame(cbind(history, 1:nrow(history)))
colnames(df_history) <- c("theta0", "theta1", "cost", "iteration")

p3 <- ggplot(df_history, aes(theta0, theta1)) + geom_path() +
  scale_x_continuous("Theta 0") + 
  scale_y_continuous("Theta 1") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position="none") +
  ggtitle("Parameters theta")
p3

#cost function
p4 <- ggplot(df_history, aes(iteration, cost)) + geom_path() +
  scale_x_continuous("Iterations") + 
  scale_y_continuous("Cost") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position="none") +
  ggtitle("Cost function")
p4



#3 - MCMC - Metropolis Hastings
set.seed(1984)
x <- df_hrv_age$age
y <- df_hrv_age$mean_hrv

#define likelihood and posterior, then multiply them to get the posterior
#likelihood
likelihood <- function(beta1, beta0, sigma){
  pred <- beta1*x + beta0
  likelihood <- sum(dnorm(y, mean = pred, sd = sigma, log = TRUE))
  return(likelihood)   
}
# Prior distribution (uniforms)
prior <- function(beta1, beta0, sigma){
  prior <- sum(dunif(beta1, min=-200, max=200, log = TRUE), 
               dunif(beta0, min=-200, max=200, log = TRUE),
               dunif(sigma, min=0, max=10000, log = TRUE))
  return (prior)
}
posterior <- function(beta1, beta0, sigma){
  return (likelihood(beta1, beta0, sigma) + prior(beta1, beta0, sigma))
}

#Metropolis-Hasting algorithm
#initialize parameters
startvalue = c(1, 1, 1)
iterations <- 100000
chain = array(dim = c(iterations+1, 3)) #3 parameters
chain[1,] = startvalue
for (i in 1:iterations){
  #draw from a normal proposal function the three elements, based on the current position
  proposal <- rnorm(3, mean = c(chain[i, 1], chain[i, 2], chain[i, 3]), sd= c(0.2, 1, 0.5))
  acceptance = exp(posterior(proposal[1], proposal[2], proposal[3]) - #posterior value at the new parameters values
                 posterior(chain[i, 1], chain[i, 2], chain[i, 3])) #posterior value at the current parameters values
  if (runif(1) < acceptance){ #on the long run, make the move with probability "acceptance"
    chain[i+1,] <- proposal #make the move
  }else{
    chain[i+1,] = chain[i,] #stay where you are
  }
}
#remove first samples
burnIn = 10000
chain <- chain[-(1:burnIn), ]
df_chain <- data.frame(chain)
colnames(df_chain) <- c("b1", "b0", "sigma")
df_chain[, "samples"] <- 1:nrow(df_chain)

mean(df_chain$b0)
mean(df_chain$b1)

#plot parameter estimates and samples from the chain
p1 <- ggplot(df_chain, aes(b1)) + geom_histogram() + 
  scale_x_continuous("Value") + 
  scale_y_continuous("Count") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  geom_vline(xintercept = mean(df_chain$b1), col = "red") +
  ggtitle("b1")
p2 <- ggplot(df_chain, aes(b0)) + geom_histogram() + 
  scale_x_continuous("Value") + 
  scale_y_continuous("Count") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  geom_vline(xintercept = mean(df_chain$b0), col = "red") +
  ggtitle("b0")
p3 <- ggplot(df_chain, aes(sigma)) + geom_histogram() + 
  scale_x_continuous("Value") + 
  scale_y_continuous("Count") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  geom_vline(xintercept = mean(df_chain$sigma), col = "red") +
  ggtitle("sigma")

p4 <- ggplot(df_chain, aes(samples, b1)) + geom_line() + 
  scale_x_continuous("Samples") + 
  scale_y_continuous("b1") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  ggtitle("Chain samples - b1")
p5 <- ggplot(df_chain, aes(samples, b0)) + geom_line() + 
  scale_x_continuous("Samples") + 
  scale_y_continuous("b0") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  ggtitle("Chain samples - b0")
p6 <- ggplot(df_chain, aes(samples, sigma)) + geom_line() + 
  scale_x_continuous("Samples") + 
  scale_y_continuous("sigma") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  ggtitle("Chain samples - sigma")

#sampler exploration of the parameter's space
ggplot(df_chain, aes(b0, b1)) + geom_path() + 
  scale_x_continuous("b0") + 
  scale_y_continuous("b1") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1))





####4 - JAGS
library(rjags)
library(coda)
set.seed(1984)

#parameters definition
x <- df_hrv_age$age
y <- df_hrv_age$mean_hrv
N <- length(y)

#estimation
jags <- jags.model(paste(dir.wd, "model_jags.R", sep = ""),
                   data = list("y" = y,
                               "x" = x,
                               "N" = N),
                   n.chains = 4,
                   n.adapt = 1000)

mcmc_samples <- coda.samples(jags,
                             c("b0", "b1", "sigma"),
                             5000)

ggplot(df_hrv_age, aes(age, mean_hrv)) + geom_point() + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Heart Rate Variability (rMSSD)") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  geom_abline(intercept = summary(mcmc_samples)$statistics[1, 1], slope = summary(mcmc_samples)$statistics[2, 1], col = "darkred") +
  ggtitle("MCMC, Gibbs sampling")












