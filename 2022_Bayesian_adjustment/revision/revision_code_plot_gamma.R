# Load the ggplot2 library for nicer graphics (optional)
library(ggplot2)
library(gridExtra)
library(grid)

# Define the shape and rate parameters for the Gamma distribution

# Create a list of vectors, where each vector contains shape and rate
params_list <- list(
  c(shape=1.1051, rate=0.1051, mode=1.0, sd=10, xrange1=0, xrange2=200,by=1),
  c(shape=1.01005, rate=0.01005, mode=1.0, sd=100,xrange1=0, xrange2=200,by=1),
  c(shape=101.99, rate=100.99, mode=1.0, sd=0.1,xrange1=0, xrange2=200,by=1)
  )



plist=list()
# Create a sequence of values from 0 to 20 (for example)

# Loop through each vector in the list
for (params in params_list) {
  shape <- params['shape']
  rate <- params['rate']
  xrange1 <- params['xrange1']
  xrange2 <- params['xrange2']
  mode <- params['mode']
  sd <- params['sd']
  by <- params['by']
  #Mode=(shape-1)/rate
  #SD=sqrt(shape)/rate
    x_values <- seq(xrange1,xrange2, by = by)
    y_values <- dgamma(x_values, shape, rate)
    gamma_data <- data.frame(x_values, y_values)
p=ggplot(gamma_data, aes(x = x_values, y = y_values)) +
  geom_line(color = "blue") +
  ggtitle(sprintf("Gamma (Mode=%s, SD=%s)", mode, sd)) +
  xlab("Value") + ylab("Density")+
  theme_minimal()
plist[[length(plist)+1]]=p
}
p=grid.arrange(grobs=plist,ncol=3)

ggsave('F:/Google Drive/post-doc/Bayesian_Project/new_model/revision/figures/gamma.png',p, 
      bg='white', width=10, height=5, units='in', dpi=300)



