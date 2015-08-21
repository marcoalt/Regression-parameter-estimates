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
plot_top <- ggplot(df_hrv_age, aes(age)) + geom_density(alpha = 0.5) + 
  scale_x_continuous("Age") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position = "none")  

# marginal density of y - plot on the right
plot_right <- ggplot(df_hrv_age, aes(mean_hrv)) + geom_density(alpha = 0.5) + 
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

#plot
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
  #scale_color_manual(values = c("orange", "purple")) + 
  scale_x_continuous("Theta 0") + 
  #stat_smooth(method = "lm") +
  scale_y_continuous("Theta 1") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position="none") +
  ggtitle("Parameters theta")
#parameters space
p3

#cost function
p4 <- ggplot(df_history, aes(iteration, cost)) + geom_path() +
  #scale_color_manual(values = c("orange", "purple")) + 
  scale_x_continuous("Iterations") + 
  #stat_smooth(method = "lm") +
  scale_y_continuous("Cost") + 
  theme(panel.background = element_rect(fill = 'ghostwhite', colour = 'ghostwhite')) +
  theme(plot.margin = unit(c(margin_plots,margin_plots,margin_plots,margin_plots), "cm")) +
  theme(legend.position="none") +
  ggtitle("Cost function")
p4




As we can see from the results obtained in this post, especially for simple toy examples, 
the three methods return basically the same solution. Especially considering the never ending 
discussions between frequentists and Bayesians, deciding which way to go might be more of a philosophical 
question than anything else. However, there are clear advantages and disadvantages to each one of 
these techniques, that I will try to summarize here:
  
  
  CHECK ALSO MY NOTES ON THE NOTEBOOK


