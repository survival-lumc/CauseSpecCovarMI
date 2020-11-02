
# Using https://www.rdatagen.net/post/a-hidden-process-part-2-of-2/

library(MASS)
library(tidyverse)

set.seed(123)
n <- 20000
Z <- rnorm(n)
beta_Z <- 1
X <- beta_Z * Z + rlogis(n) # standard logistic error
# Same as X <- rlogis(n, location = beta_Z * Z)

# Try with 2 categories

breaks_X <- c(-Inf, -1, 1.5, Inf)
catX <- cut(X, breaks = breaks_X, ordered_result = T)
breakos <- data.frame(breaks_X, "dens" = dlogis(breaks_X))

dat <- data.frame(catX, Z)

?polr
mod <- MASS::polr(catX ~ Z, data = dat)
summary(mod)
cor(X, Z)

mod$zeta # log odds of being in that cat or lower when Z = 0
mod$coefficients # Param as -beta, same across classes 

predict(mod, newdata = data.frame(Z = 0), type = "p")
plogis(mod$zeta)
diff(c(0, plogis(mod$zeta), 1))

predict(mod, newdata = data.frame(Z = 1), type = "p")
diff(c(0, plogis(mod$zeta - mod$coefficients * 1), 1))



dat %>% 
  ggplot(aes(Z, X)) +
  geom_point() +
  geom_hline(data = breakos[2:3, ], aes(yintercept = breaks_X),
             col = "blue", linetype = "dashed")

Z_test <- 0:3


lapply(Z_test, function(i) {
  X0 <- X + i * beta_Z
  cbind.data.frame(
    "Z" = i,
    "X" = X0,
    "cato" = cut(X0, breaks = breaks_X, ordered_result = T),
    "dens" = dlogis(X0, location = i * beta_Z)
  )
}) %>% 
  bind_rows() %>% 
  mutate(Z = as.factor(Z)) %>% 
  ggplot(aes(X, dens)) +
  geom_line() +
  coord_flip() +
  facet_wrap(~ Z, nrow = 1) +
  geom_area(aes(fill = cato, group = cato))
  



dat %>% 
  mutate(dens = dlogis(X)) %>% 
  ggplot(aes(X, dens)) +
  geom_line() +
  geom_segment(data = breakos,
               aes(x = breaks_X, xend = breaks_X,
                   y = 0, yend = dens), linetype = "dashed")


dat %>% 
  mutate(dens = dlogis(X)) %>% 
  ggplot(aes(X, dens)) +
  geom_line() +
  geom_segment(data = breakos,
               aes(x = breaks_X - beta_Z * 1, xend = breaks_X - beta_Z * 1,
                   y = 0, yend = dlogis(breaks_X - beta_Z * 1)), linetype = "dashed")






eta <- model.matrix(~ factor(catX, ordered = F) + Z) %*% c(0, 0.1, 0.3, 1)

# Generate survival data here...
t <- rexp(n, rate = 0.2 * exp(eta))
hist(t)

library(survival)
data <- cbind.data.frame("catX" = factor(catX, ordered = T), X, Z, t, "eps" = 1)

options("contrasts")
mod1 <- coxph(Surv(t, eps) ~ catX + Z, data = data)
mod1

coxph(Surv(t, eps) ~ C(catX, contr = "contr.treatment") + Z, data = data)
coxph(Surv(t, eps) ~ factor(catX, ordered = F) + Z, data = data)

miss_ind <- rbinom(n, 1, prob = plogis(data$Z))
data$X <- ifelse(miss_ind == 1, NA, data$catX)
data$X <- factor(data$X, levels = c(1, 2, 3),
                 labels = c("Low", "Med", "High"), ordered = T)
data$miss_ind = miss_ind

data %>% 
  ggplot(aes(Z, catX, col = factor(miss_ind))) +
  geom_jitter()

library(survival)
library(smcfcs)
mod <- smcfcs(
  originaldata = data,
  smtype = "coxph", 
  smformula = "Surv(t, eps) ~ C(X, contr = 'contr.treatment') + Z",
  m = 5,
  method = c("", "", "", "", "podds", "")
)

# General defaults - this will solve all harms
options(contrasts = rep("contr.treatment", 2))  

mod <- smcfcs(
  originaldata = data,
  smtype = "coxph", 
  smformula = "Surv(t, eps) ~ X + Z",
  m = 10, numit = 25,
  method = c("", "", "", "", "podds")
)

mod <- smcfcs(
  originaldata = data,
  smtype = "coxph", 
  smformula = "Surv(t, eps) ~ X + Z",
  m = 5,
  method = c("", "", "", "", "podds")
)


purrr::map_dfr(1:10, function(i) {
  as.data.frame(t(mod$smCoefIter[i, ,])) %>% 
    mutate(iter = 1:n())
}, .id = "imp_dat") %>% 
  rename("beta_med" = V1, "beta_high" = V2,
         "beta_Z" = V3) %>% 
  gather("var", "value", beta_med:beta_Z) %>% 
  ggplot(aes(iter, value, col = imp_dat)) +
  #stat_summary(aes(group = 1), fun = mean,
  #             geom = "line", size = 1.5, col = "black") + 
  geom_line() +
  theme(legend.position = "none") +
  ggplot2::facet_wrap(~ var) 
  





mod2.1 <- lm(Z ~ X, data = data)
mod2.2<- lm(Z ~ C(X, contr = 'contr.treatment'), data = data)
mod2.1
mod2.2

predict(mod2.2)
predict(mod2.2,  newdata = data.frame(X = c("Low", "Med", "High", NA)))
#, "Med", "High", NA)))
predict(mod2.2,  newdata = data.frame(X = factor(c("Low", "Med", "High", NA),
                                                 levels = c("Low", "Med", "High" NA),
                                                 ordered = T)))



predict(mod2.1, newdata = data.frame(X = "Med")) 
predict(mod2.2,  newdata = data.frame(X = c("Low", "Med")))
predict(mod2.2,  newdata = data.frame(X = c("Low", "Med", "High", NA)))


# Great sources of contrasts:
# https://stats.stackexchange.com/questions/105115/polynomial-contrasts-for-regression
# https://link.springer.com/chapter/10.1007/978-1-4612-0971-3_10


# Try this same code for logistic
# https://www.rdatagen.net/post/ordinal-regression/
        
table(catX)




# Demo 01/07/2020 ----------------------------------------------------------

head(data)

dat %>% 
  ggplot(aes(Z, X)) +
  geom_point(col = "darkblue") +
  geom_hline(data = breakos[2:3, ], aes(yintercept = breaks_X),
             col = "black", linetype = "dashed", size = 1.5)

data %>% 
  ggplot(aes(catX, Z)) +
  geom_jitter(width = 0.2, col = "darkblue", alpha = 0.5)

data$catX

coxph(Surv(t, eps) ~ catX + Z, data = data)

options("contrasts")
contr.treatment(n = 3)
contr.poly(n = 3)

coxph(Surv(t, eps) ~ C(catX, contr = "contr.treatment") + Z, data = data)




options(contrasts = rep("contr.treatment", 2))  
options("contrasts")

coxph(Surv(t, eps) ~ catX + Z, data = data)


miss_ind <- rbinom(n, 1, prob = plogis(data$Z))
data$X <- ifelse(miss_ind == 1, NA, data$catX)
data$X <- factor(data$X, levels = c(1, 2, 3),
                 labels = c("Low", "Med", "High"), ordered = T)
data$miss_ind = miss_ind

data %>% 
  ggplot(aes(Z, catX, col = factor(miss_ind))) +
  geom_jitter()

mod <- smcfcs(
  originaldata = data,
  smtype = "coxph", 
  smformula = "Surv(t, eps) ~ X + Z",
  m = 5,
  method = c("", "", "", "", "podds")
)
