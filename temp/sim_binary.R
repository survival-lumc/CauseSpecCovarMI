library(MASS)

# pick pick proportion first, then point of dichot
p <- .25
x0 <- qnorm(p, lower.tail = F)

# https://www.tandfonline.com/doi/full/10.1080/03610918.2012.707725
# Section 2.3
# or also # https://stats.stackexchange.com/questions/313861/generate-a-gaussian-and-a-binary-random-variables-with-predefined-correlation
# This will work for standard normal data
pbiserial_to_pearson <- function(p, x0, r_pb) {
  
  h <- dnorm(x0) # 'ordinate' = density/height of curve
  r <- r_pb * sqrt(p * (1 - p)) / h
  return(r)
}


r <- pbis_to_pears(p, x0, r_pb = .5)


dat <- mvrnorm(n = 1000, mu = c(0, 0), Sigma = matrix(c(1, r, 
                                                        r, 1), nrow = 2))


dat <- as.data.frame(dat)
colnames(dat) <- c("X", "Z")
cor(dat$X, dat$Z)
dat$X_bin <- with(dat, ifelse(X <= x0, 0, 1))
cor(dat$X_bin, dat$Z)

#plot(dat$X, dat$Z)
#plot(dat$X_bin, dat$Z)

lmod <- with(dat, lm(X_bin ~ Z))

logreg <- with(dat, glm(X_bin ~ Z, family = binomial))

plot(logreg)

table(dat$X_bin)
# Set 1 - p = 0.5
x0 <- qnorm(0.5)

rho <- 0.5

r <- (rho * sqrt(0.5^2)) / dnorm(x0)
r





