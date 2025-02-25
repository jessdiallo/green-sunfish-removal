# Supplementary material for Diallo et al. manuscript:

# Diallo, Jessica O., Sarah J. Converse, Matthew Chmiel, Andy Stites, & Julian D. Olden. 2025. "Optimizing control of a freshwater invader in time and space." Ecological Applications.

# Corresponding author: Jessica O. Diallo
# Email: jodiallo@uw.edu; jessicaodiallo@gmail.com

# load packages
library(tidyverse)
library(rjags)
library(jagsUI)
library(cowplot)
library(scales)

set.seed(30)

# load data
load("SM-data.RData")

# JAGS Model - McGee Wash ----
sink("model.jags.MGW")
cat("
model {

for(i in 1:I) {
  for(k in 1:K) {
    p.miss[i,k,1:700] <- rep(1/700,700)
    miss[i,k] ~ dcat(p.miss[i,k,1:700])
    N[i,1,k,1] <-  N.base[i,1,k,1] + miss[i,k]
  }
}

beta.p ~ dunif(0, 10)
alpha.p[1] ~ dunif(-10,10)
alpha.p[2] ~ dunif(-10,10)
alpha.p[3] ~ dunif(-10,10)
alpha.p[4] <- alpha.p[1]
alpha.p[5] <- alpha.p[2]
alpha.p[6] <- alpha.p[3]
s1 ~ dunif(0, 1)
s2 ~ dunif(0, 1)

for(i in 1:4){
  f[i] ~ dnorm(0,1/(12^2)) T(0,) 
}

for(t in 1:T) {
  for(i in 1:I){
    for(j in 1:J){ 
      logit(p[i,j,t]) <- alpha.p[j] + beta.p * log(xi[i,j,t])
    }
    for(k in 1:K) {
      for(j in 1:J) {                
        Y[i,j,k,t] ~ dbin(p[i,j,t], N[i,j,k,t])
      }
      for(j in 2:J) {
        N[i,j,k,t] <- N[i,j-1,k,t] - Y[i,j-1,k,t]
      }
      R[i,k,t] <- N[i,J,k,t] - Y[i,J,k,t]
    }
  }
}

for(t in 1:(T-1)){
  P[1,1,t] <- ifelse(t == 11 || t == 23 || t == 33 || t == 41, 0,s1^days_btwn[t])
  P[2,1,t] <- ifelse(t == 11 || t == 23 || t == 33 || t == 41, s1^days_btwn[t],0)
  P[2,2,t] <- s2^days_btwn[t]
}

for(t in 1:12){
  P[1,2,t] <- ifelse(t == 11, f1,0)
}
for(t in 13:23){
  P[1,2,t] <- ifelse(t == 23, f2,0)
}
for(t in 24:33){
  P[1,2,t] <- ifelse(t == 33, f3,0)
}
for(t in 34:(T-1)){
  P[1,2,t] <- ifelse(t == 41, f4,0)
}

for(t in 1:(T-1)){
  for(i in 1:I){
    D[i,1,t] ~ dpois(R[i,1,t]*P[1,1,t] + R[i,2,t]*P[1,2,t]) # juveniles
    D[i,2,t] ~ dpois(R[i,1,t]*P[2,1,t] + R[i,2,t]*P[2,2,t]) # adults
    for(k in 1:K){
      N[i,1,k,t+1] <- D[i,k,t]
    }
  }
}

for(t in 1:T){
  N.sum[t] <- sum(N[,1,,t])
}

} # end model
", fill= TRUE)
sink()

I <- 26 # number of pools
J <- 6 # number of gear types
K <- 2 # number of size classes
T <- 42 # time period

# Initial values
Y <- array_count_MGW
N.base <- array(NA_real_, dim = c(I, J, K, T))
for(i in 1:I){
  N.base[i,1,1,1] <- sum(Y[i,,1,1:10])
  N.base[i,1,2,1] <- sum(Y[i,,2,1:10])
}

# initialize D to > than the number that will be removed in the following year
Y <- array_count_MGW
Y.remove1 <- apply(Y[,,1,],c(1,3),sum)
Y.remove2 <- apply(Y[,,2,],c(1,3),sum)

D.init <- array(NA,dim = c(I,K,(T-1)))
for(t in 1:(T-1)){
  D.init[,1,t] <- Y.remove1[,t+1] + 1
  D.init[,2,t] <- Y.remove2[,t+1] + 1
}

f.init <- c(12,0.5,12,0.5)

inits <- function (){
  list(D=D.init, f=f.init)
}

# Bundle data together
data <- list(Y = array_count_MGW,
             I = 26, J = 6, K = 2, T = 42, 
             xi = array_stdeffort_MGW, 
             days_btwn = date_diff_MGW, N.base = N.base)

# Parameters monitored
parameters <- c("N", "N.sum", "alpha.p", "beta.p", "s1", "s2", 
                "f1", "f2", "f3", "f4", "p")

# MCMC settings
ni <- 20000
nt <- 1
nb <- 10000
nc <- 3

# Call JAGS Function
output_MGW <- jags(data, 
                   inits, 
                   parameters, 
                   "model.jags.MGW", 
                   n.chains = nc, 
                   n.thin = nt, 
                   n.iter = ni,
                   n.burnin = nb)


# JAGS Model - Eash Ash Creek ----
sink("model.jags.EAC")
cat("
model {

for(i in 1:I) {
  for(k in 1:K) {
    p.miss[i,k,1:700] <- rep(1/700,700)
    miss[i,k] ~ dcat(p.miss[i,k,1:700])
    N[i,1,k,1] <-  N.base[i,1,k,1] + miss[i,k]
  }
}
beta.p ~ dunif(0, 10)
alpha.p[1] ~ dnorm(alpha2mean, 1/(alpha2sd^2))
alpha.p[2] ~ dnorm(alpha3mean, 1/(alpha3sd^2))
alpha.p[3] <- alpha.p[1]
alpha.p[4] <- alpha.p[2]
s1 ~ dnorm(s1mean, 1/(s1sd^2))      # priors for survival (daily)
s2 ~ dnorm(s2mean, 1/(s2sd^2))

f ~ dnorm(0,1/(12^2)) T(0,)

for(t in 1:T) {
  for(i in 1:I){
    for(j in 1:J){ 
      logit(p[i,j,t]) <- alpha.p[j] + beta.p * log(xi[i,j,t])
    }
    for(k in 1:K) {
      for(j in 1:J) {                
        Y[i,j,k,t] ~ dbin(p[i,j,t], N[i,j,k,t])
      }
      for(j in 2:J) {
        N[i,j,k,t] <- N[i,j-1,k,t] - Y[i,j-1,k,t]
      }
      R[i,k,t] <- N[i,J,k,t] - Y[i,J,k,t]
    }  
  }
}

for(t in 1:(T-1)){
  P[1,2,t] <- ifelse(t == 8, f,0)
  P[1,1,t] <- ifelse(t == 8, 0,s1^days_btwn[t])
  P[2,1,t] <- ifelse(t == 8, s1^days_btwn[t],0)
  P[2,2,t] <- s2^days_btwn[t]
  for(i in 1:I){
    for(k in 1:K){
      D[i,k,t] ~ dpois(R[i,1,t]*P[k,1,t] + R[i,2,t]*P[k,2,t])
    }
    for(k in 1:K){
      N[i,1,k,t+1] <- D[i,k,t] 
    }
  }
}

for(t in 1:T){
  N.sum[t] <- sum(N[,1,,t])
}

} # end model
", fill= TRUE)
sink()

I <- 39 # number of pools
J <- 4 # number of gear types
K <- 2 # number of size classes
T <- 17 # time period

# Initial values
Y <- array_count_EAC
N.base <- array(NA_real_, dim = c(I, J, K, T))
for(i in 1:I){
  N.base[i,1,1,1] <- sum(Y[i,,1,1:7])
  N.base[i,1,2,1] <- sum(Y[i,,2,1:7])
}

# initialize D to >> than the number that will be removed in the following year
Y.remove1 <- apply(Y[,,1,],c(1,3),sum)
Y.remove2 <- apply(Y[,,2,],c(1,3),sum)

D.init <- array(NA,dim = c(I,K,(T-1)))
for(t in 1:(T-1)){
  D.init[,1,t] <- Y.remove1[,t+1] + 1
  D.init[,2,t] <- Y.remove2[,t+1] + 1
}

inits <- function (){
  list(D=D.init)
}

# Bundle data together
data <- list(Y = array_count_EAC,
             I = 39, J = 4, K = 2, T = 17,
             xi = array_stdeffort_EAC, 
             days_btwn = date_diff_EAC, N.base = N.base,
             alpha2mean = MGW_parameters$alpha.2[1],
             alpha2sd = MGW_parameters$alpha.2[2],
             alpha3mean = MGW_parameters$alpha.3[1],
             alpha3sd = MGW_parameters$alpha.3[2],
             s1mean = MGW_parameters$s1[1],
             s1sd = MGW_parameters$s1[2],
             s2mean = MGW_parameters$s2[1],
             s2sd = MGW_parameters$s2[2])

# Parameters monitored
parameters <- c("N", "N.sum", "alpha.p", "beta.p", "s1", "s2", "f", "p")

# MCMC settings
ni <- 20000
nt <- 1
nb <- 10000
nc <- 3


# Call JAGS Function
output_EAC <- jags(data, 
                   inits, 
                   parameters, 
                   "model.jags.EAC", 
                   n.chains = nc, 
                   n.thin = nt, 
                   n.iter = ni,
                   n.burnin = nb)

# Simulations ----
## Spatial, Random (Fig 5) ----
# Random Selection
I <- 26 # number of pools
J <- 6 # number of passes/gear types
K <- 2 # number of size classes
T <- 42 # removal number (years)
Sim <- 1000 # number of simulations
M <- 11

# random selection of 1000 iterations from 30,000
param <- rdunif(1000,1,30000)

# Random selection set-up
percentage <- seq(0,100,10)
selection <- round(percentage/100*I)

# Starting arrays and values
Y <- array(0, dim = c(I,J,K,T,Sim,M))
N <- array(0, dim = c(I,J,K,T,Sim,M))
R <- array(0, dim = c(I,K,T,Sim,M))
D <- array(0, dim = c(I,K,T,Sim,M))
L <- array(0, dim = c(K,K,T,Sim)) # transition matrix
p <- output_p_MGW # probability of capture dim = c(26,6,42)

# Dates
# original object: dates_MGW
days_btwn <- date_diff_MGW # days between removals
T_hatch <- c(11, 23, 32, 41)

for(s in 1:Sim){
  for(t in 1:T) {
    if(t %in% T_hatch){  # these time points allowing new size 1 individuals 
      L[1,1,t,s] <- 0
      L[1,2,t,s] <- f1.fullarray_MGW$f1[param[s]]
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
    } else {
      L[1,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
      L[1,2,t,s] <- 0
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- 0
    }
  }
}

for(m in 1:M){
  for(s in 1:Sim){
    for(i in 1:I){
      for(k in 1:K){
        N[i,1,k,1,s,m] <- starting_N[i,k]
      }
    }
    for(t in 1:T) {
      removal_pools <- sample(1:I, selection[m])
      for(i in 1:I){
        for(k in 1:K){
          for(j in 1:J){
            if(i %in% removal_pools){
              Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
            } else {
              Y[i,j,k,t,s,m] <- 0
            }
            if(j < J) {
              N[i,j+1,k,t,s,m] <- N[i,j,k,t,s,m] - Y[i,j,k,t,s,m]
            }
          } # end of J for loop
          R[i,k,t,s,m] <- N[i,J,k,t,s,m] - Y[i,J,k,t,s,m] 
        } # end of K for loop 
        if(t<T) { # Not including t = T here because with last time step there isn't
          # population growth and movement
          for(k in 1:K){
            # Population growth
            D[i,k,t,s,m] <- trunc(R[i,1,t,s,m]*L[k,1,t,s] + R[i,2,t,s,m]*L[k,2,t,s])
            N[i,1,k,t+1,s,m] <- D[i,k,t,s,m]
          } # end of k for loop
        } # end of if t<T
      } # end of I for loop
    } # end of T for loop
  } # end of Sim for loop
} # end of M for loop


# set up data frame with simulation output
# sum matrices
N_array <- array(NA, dim = c(T,Sim,M))
Y_array <- array(NA, dim = c(T,Sim,M))

for(m in 1:M){
  for(s in 1:Sim){
    for(t in 1:T){
      N_array[t,s,m] <- sum(N[,1,,t,s,m])
      Y_array[t,s,m] <- sum(Y[,1:3,,t,s,m])
    }
  }
}

# Just fish at T=42
N_matrix_t42 <- N_array[42,,]

M_col2 <- vector()
for(m in 1:M){
  M_col2 <- c(M_col2, rep(m, times = Sim))
}

N_t42 <- vector()
for(m in 1:M){
  for(s in 1:Sim){
    N_t42 <- c(N_t42, N_matrix_t42[s,m])
  }
}

N_t42_df1 <- tibble(M_col2, N_t42)
N_t42_df1$M_col2 <- factor(N_t42_df1$M_col2, levels = seq(1,11,1))

# Random, Continuous Selection
I <- 26 # number of pools
J <- 6 # number of passes/gear types
K <- 2 # number of size classes
T <- 42 # removal number (years)
Sim <- 1000 # number of simulations
M <- 11

# Random selection set-up
percentage <- seq(0,100,10)
selection <- round(percentage/100*I)

# Starting arrays and values
Y <- array(0, dim = c(I,J,K,T,Sim,M))
N <- array(0, dim = c(I,J,K,T,Sim,M))
R <- array(0, dim = c(I,K,T,Sim,M))
D <- array(0, dim = c(I,K,T,Sim,M))
L <- array(0, dim = c(K,K,T,Sim)) # transition matrix
p <- output_p_MGW # probability of capture dim = c(26,6,42)

# random selection of 1000 iterations from 30,000
param <- rdunif(1000,1,30000)

# Dates
# original object: dates_MGW
days_btwn <- date_diff_MGW # days between removals
T_hatch <- c(11, 23, 32, 41)

for(s in 1:Sim){
  for(t in 1:T) {
    if(t %in% T_hatch){  # these time points allowing new size 1 individuals 
      L[1,1,t,s] <- 0
      L[1,2,t,s] <- f1.fullarray_MGW$f1[param[s]]
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
    } else {
      L[1,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
      L[1,2,t,s] <- 0
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- 0
    }
  }
}

for(m in 1:M){
  for(s in 1:Sim){
    for(i in 1:I){
      for(k in 1:K){
        N[i,1,k,1,s,m] <- starting_N[i,k]
      }
    }
    for(t in 1:T) {
      if(selection[m] > 0){
        first_removal_pool <- sample(1:(I-(selection[m]-1)), 1)
        removal_pools <- seq(from = first_removal_pool, 
                             to = (first_removal_pool + (selection[m]-1)))
      } else {
        removal_pools <- 0
      }
      for(i in 1:I){
        for(k in 1:K){
          for(j in 1:J){
            if(i %in% removal_pools){
              Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
            } else {
              Y[i,j,k,t,s,m] <- 0
            }
            if(j < J) {
              N[i,j+1,k,t,s,m] <- N[i,j,k,t,s,m] - Y[i,j,k,t,s,m]
            }
          } # end of J for loop
          R[i,k,t,s,m] <- N[i,J,k,t,s,m] - Y[i,J,k,t,s,m] 
        } # end of K for loop 
        if(t<T) { # Not including t = T here because with last time step there isn't
          # population growth and movement
          for(k in 1:K){
            # Population growth
            D[i,k,t,s,m] <- trunc(R[i,1,t,s,m]*L[k,1,t,s] + R[i,2,t,s,m]*L[k,2,t,s])
            N[i,1,k,t+1,s,m] <- D[i,k,t,s,m]
          } # end of k for loop
        } # end of if t<T
      } # end of I for loop
    } # end of T for loop
  } # end of Sim for loop
} # end of P for loop

# set up data frame with simulation output
# sum matrices
N_array <- array(NA, dim = c(T,Sim,M))
Y_array <- array(NA, dim = c(T,Sim,M))

for(m in 1:M){
  for(s in 1:Sim){
    for(t in 1:T){
      N_array[t,s,m] <- sum(N[,1,,t,s,m])
      Y_array[t,s,m] <- sum(Y[,1:3,,t,s,m])
    }
  }
}

# Just fish at T=42
# use N_array, not N_summary b/c I need all the simulations for the boxplots
N_matrix_t42 <- N_array[42,,]

M_col2 <- vector()
for(m in 1:M){
  M_col2 <- c(M_col2, rep(m, times = Sim))
}

N_t42 <- vector()
for(m in 1:M){
  for(s in 1:Sim){
    N_t42 <- c(N_t42, N_matrix_t42[s,m])
  }
}

N_t42_df2 <- tibble(M_col2, N_t42)
N_t42_df2$M_col2 <- factor(N_t42_df2$M_col2, levels = seq(1,11,1))

# Saving first two spatial together
# N_t42_df1 and N_t42_df2
N_t42_df1$sim <- "Random"
N_t42_df2$sim <- "Random Contiguous"
# sim_output_spatial_random <- bind_rows(N_t42_df1, N_t42_df2)

## Spatial, Volume (Fig A2) ----
I <- 26 # number of pools
J <- 6 # number of passes/gear types
K <- 2 # number of size classes
T <- 42 # removal number (years)
Sim <- 1000 # number of simulations
M <- 11

# random selection of 1000 iterations from 30,000
param <- rdunif(1000,1,30000)

# Random selection set-up
percentage <- seq(0,100,10)
selection <- round(percentage/100*I)

# Starting arrays and values
Y <- array(0, dim = c(I,J,K,T,Sim,M))
N <- array(0, dim = c(I,J,K,T,Sim,M))
R <- array(0, dim = c(I,K,T,Sim,M))
D <- array(0, dim = c(I,K,T,Sim,M))
L <- array(0, dim = c(K,K,T,Sim)) # transition matrix
p <- output_p_MGW # probability of capture dim = c(26,6,42)

# Dates
# original object: dates_MGW
days_btwn <- date_diff_MGW # days between removals
T_hatch <- c(11, 23, 32, 41)

# pool volume
colnames(pool_data) <- c("site", "easting", "northing", "photo_ID", "max_depth",
                         "length", "width", "tail_crest_depth", "head_crest_depth")
pool_dim <- pool_data %>%
  select(site, max_depth, length, width)
volume <- (1/2)*(4/3)*(3.1415926535) * pool_dim$max_depth * (1/2)*pool_dim$length * (1/2)*pool_dim$width
pool_vol <- cbind(pool_dim, volume)
pool_vol <- pool_vol %>%
  select(site, volume)
pool_vol <- as_tibble(pool_vol) %>%
  arrange(-volume)
pool_vol$site2 <- as.numeric(substr(pool_vol$site, 6,7))

for(s in 1:Sim){
  for(t in 1:T) {
    if(t %in% T_hatch){  # these time points allowing new size 1 individuals 
      L[1,1,t,s] <- 0
      L[1,2,t,s] <- f1.fullarray_MGW$f1[param[s]]
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
    } else {
      L[1,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
      L[1,2,t,s] <- 0
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- 0
    }
  }
}

for(m in 1:M){
  for(s in 1:Sim){
    for(i in 1:I){
      for(k in 1:K){
        N[i,1,k,1,s,m] <- starting_N[i,k]
      }
    }
    for(t in 1:T) {
      if(selection[m] > 0){
        removal_pools <- pool_vol$site2[1:selection[m]]
      } else {
        removal_pools <- 0
      }
      for(i in 1:I){
        for(k in 1:K){
          for(j in 1:J){
            if(i %in% removal_pools){
              Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
            } else {
              Y[i,j,k,t,s,m] <- 0
            }
            if(j < J) {
              N[i,j+1,k,t,s,m] <- N[i,j,k,t,s,m] - Y[i,j,k,t,s,m]
            }
          } # end of J for loop
          R[i,k,t,s,m] <- N[i,J,k,t,s,m] - Y[i,J,k,t,s,m] 
        } # end of K for loop 
        if(t<T) { # Not including t = T here because with last time step there isn't
          # population growth and movement
          for(k in 1:K){
            # Population growth
            D[i,k,t,s,m] <- trunc(R[i,1,t,s,m]*L[k,1,t,s] + R[i,2,t,s,m]*L[k,2,t,s])
            N[i,1,k,t+1,s,m] <- D[i,k,t,s,m]
          } # end of k for loop
        } # end of if t<T
      } # end of I for loop
    } # end of T for loop
  } # end of Sim for loop
} # end of P for loop


# set up data frame with simulation output
# sum matrices
N_array <- array(NA, dim = c(T,Sim,M))
Y_array <- array(NA, dim = c(T,Sim,M))

for(m in 1:M){
  for(s in 1:Sim){
    for(t in 1:T){
      N_array[t,s,m] <- sum(N[,1,,t,s,m])
      Y_array[t,s,m] <- sum(Y[,1:3,,t,s,m])
    }
  }
}

# Just fish at T=42
# use N_array, not N_summary b/c I need all the simulations for the boxplots
N_matrix_t42 <- N_array[42,,]

M_col2 <- vector()
for(m in 1:M){
  M_col2 <- c(M_col2, rep(m, times = Sim))
}

N_t42 <- vector()
for(m in 1:M){
  for(s in 1:Sim){
    N_t42 <- c(N_t42, N_matrix_t42[s,m])
  }
}

N_t42_df3 <- tibble(M_col2, N_t42)
N_t42_df3$M_col2 <- factor(N_t42_df3$M_col2, levels = seq(1,11,1))
# sim_output_spatial_vol <- N_t42_df3

## Temporal, Fewer (Fig 6) ----
# Skipping dates
I <- 26 # number of pools
J <- 6 # number of passes/gear types
K <- 2 # number of size classes
T <- 42 # removal number (years)
Sim <- 1000 # number of simulations
M <- 8

# random selection of 1000 iterations from 30,000
param <- rdunif(1000,1,30000)

# Starting arrays and values
Y <- array(0, dim = c(I,J,K,T,Sim,M))
N <- array(0, dim = c(I,J,K,T,Sim,M))
R <- array(0, dim = c(I,K,T,Sim,M))
D <- array(0, dim = c(I,K,T,Sim,M))
L <- array(0, dim = c(K,K,T,Sim)) # transition matrix
p <- output_p_MGW # probability of capture dim = c(26,6,42)

# Dates
# original object: dates_MGW
days_btwn <- date_diff_MGW # days between removals
# post-june hatch
T_post_hatch <- c(12, 24, 34, 42)
# june hatch
T_hatch <- c(11, 23, 33, 41)
# may spawn
T_spawn <- c(10, 22, 32, 40)
# pre-may spawn
T_pre_spawn <- c(9, 21, 31, 39)

removal_index_vec1 <- c(1,12,24,34,42) # 12% (5) skipping 10,11,9,7 (Aug, June, July, July, July)
removal_index_vec2 <- c(1,6,10,15,19,24,28,34,37,42) # 23% (10) skipping 4,3,4,3,4,3,5,2,4
removal_index_vec3 <- c(1,4,7,10,13,16,19,22,24,27,30,34,37,40,42) # 36% (15) skipping 2,2,2,..3..,1
removal_index_vec4 <- c(1,3,5,7,9,12,14,16,18,20,22,24,27,29,31,34,36,38,40,42) 
# 47% (20) 2,2,2,2,3,2,2,2,2,2,2,3,2,2,3,2,2,2,2
removal_index_vec5 <- c(1,3,5,6,8,10,12,13,15,17,18,20,22,24,25,27,29,30,32,34,36,37,39,41,42) 
# 60% (25) 2,2,1,2,2,2,1,2,2,1,2,2,2,1,2,2,1,2,2,2,1,2,2,1
removal_index_vec6 <- c(1,2,4,5,6,8,9,11,12,13,15,16,18,19,20,22,23,25,26,27,29,30,32,33,34,36,37,39,40,42) 
# 71% (30) 1,2,1,1,2,1,2,1,1,2,1,2,1,1,2,1,2,1,1,2,1,2,1,1,2,1,2,1,2
removal_index_vec7 <- c(1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23,25,26,27,28,29,31,32,33,34,35,37,38,39,40,42) 
# 83% (35) 1,1,1,1,2,1,1,1,1,2,1,1,1,1,2,1,1,1,1,2,1,1,1,1,2,1,1,1,1,2,1,1,1,2
removal_index_vec8 <- seq(1,42) # 42


# create list
removal_index_list <- list(removal_index_vec1,removal_index_vec2, removal_index_vec3,
                           removal_index_vec4,removal_index_vec5, removal_index_vec6,
                           removal_index_vec7, removal_index_vec8)

for(s in 1:Sim){
  for(t in 1:T) {
    if(t %in% T_hatch){  # these time points allowing new size 1 individuals 
      L[1,1,t,s] <- 0
      L[1,2,t,s] <- f1.fullarray_MGW$f1[param[s]]
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
    } else {
      L[1,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
      L[1,2,t,s] <- 0
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- 0
    }
  }
}

for(m in 1:M){
  for(s in 1:Sim){
    for(i in 1:I){
      for(k in 1:K){
        N[i,1,k,1,s,m] <- starting_N[i,k]
      }
    }
    for(t in 1:T) {
      for(i in 1:I){
        for(k in 1:K){
          for(j in 1:J){
            if(t %in% removal_index_list[[m]]){
              Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
            } else {
              Y[i,j,k,t,s,m] <- 0
            }
            if(j < J) {
              N[i,j+1,k,t,s,m] <- N[i,j,k,t,s,m] - Y[i,j,k,t,s,m]
            }
          } # end of J for loop
          R[i,k,t,s,m] <- N[i,J,k,t,s,m] - Y[i,J,k,t,s,m] 
        } # end of K for loop
        if(t<T) { # Not including t = T here because with last time step there isn't
          # population growth and movement
          for(k in 1:K){
            # Population growth
            D[i,k,t,s,m] <- trunc(R[i,1,t,s,m]*L[k,1,t,s] + R[i,2,t,s,m]*L[k,2,t,s])
            N[i,1,k,t+1,s,m] <- D[i,k,t,s,m]
          } # end of k for loop
        } # end of if t<T
      } # end of I for loop
    } # end of T for loop
  } # end of Sim for loop
} # end of M for loop


# set up data frame with simulation output
# sum matrices
N_array <- array(NA, dim = c(T,Sim,M))
Y_array <- array(NA, dim = c(T,Sim,M))

for(m in 1:M){
  for(s in 1:Sim){
    for(t in 1:T){
      N_array[t,s,m] <- sum(N[,1,,t,s,m])
      Y_array[t,s,m] <- sum(Y[,1:3,,t,s,m])
    }
  }
}

# Just fish at T=42
# use N_array, not N_summary b/c I need all the simulations for the boxplots
N_matrix_t42 <- N_array[42,,]

M_col2 <- vector()
for(m in 1:M){
  M_col2 <- c(M_col2, rep(m, times = Sim))
}

N_t42 <- vector()
for(m in 1:M){
  for(s in 1:Sim){
    N_t42 <- c(N_t42, N_matrix_t42[s,m])
  }
}

N_t42_df4 <- tibble(M_col2, N_t42)
N_t42_df4$M_col2 <- factor(N_t42_df4$M_col2, levels = seq(1,8,1))

# Skipping dates + pre-spawn 
I <- 26 # number of pools
J <- 6 # number of passes/gear types
K <- 2 # number of size classes
T <- 42 # removal number (years)
Sim <- 1000 # number of simulations
M <- 8

# random selection of 1000 iterations from 30,000
param <- rdunif(1000,1,30000)

# Starting arrays and values
Y <- array(0, dim = c(I,J,K,T,Sim,M))
N <- array(0, dim = c(I,J,K,T,Sim,M))
R <- array(0, dim = c(I,K,T,Sim,M))
D <- array(0, dim = c(I,K,T,Sim,M))
L <- array(0, dim = c(K,K,T,Sim)) # transition matrix
p <- output_p_MGW # probability of capture dim = c(26,6,42)

# Dates
# original object: dates_MGW
days_btwn <- date_diff_MGW # days between removals
# post-june hatch
T_post_hatch <- c(12, 24, 34, 42)
# june hatch
T_hatch <- c(11, 23, 33, 41)
# may spawn
T_spawn <- c(10, 22, 32, 40)
# pre-may spawn
T_pre_spawn <- c(9, 21, 31, 39)

removal_index_vec1 <- c(1,9,21,31,39) # 12% (5) skipping 7,11,9,7
removal_index_vec2 <- c(1,5,9,14,19,21,27,31,35,39) # 23% (10) skipping 3,3,4,4,2,5,3,3,3
removal_index_vec3 <- c(1,4,7,9,13,16,19,21,24,27,31,33,37,39,42) # 36% (15) skipping 2,2,1,3,2,2,1,2,2,2,1,3,1,2
removal_index_vec4 <- c(1,3,5,7,9,10,14,16,18,21,23,25,27,29,31,33,35,37,39,42) # 47% (20) 
removal_index_vec5 <- c(1,3,5,6,8,10,12,13,15,17,18,20,22,24,25,27,29,30,31,34,36,37,39,40,42) # 60% (25) 
removal_index_vec6 <- c(1,2,4,5,6,8,9,10,12,13,15,16,18,19,20,21,23,25,26,27,29,30,31,33,34,36,37,39,40,42) # 71% (30) 
removal_index_vec7 <- c(1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23,25,26,27,28,29,31,32,33,34,35,37,38,39,40,42) # 83% (35) 
removal_index_vec8 <- seq(1,42) # 42


# create list
removal_index_list <- list(removal_index_vec1, removal_index_vec2, removal_index_vec3,
                           removal_index_vec4, removal_index_vec5, removal_index_vec6,
                           removal_index_vec7, removal_index_vec8)

for(s in 1:Sim){
  for(t in 1:T) {
    if(t %in% T_hatch){  # these time points allowing new size 1 individuals 
      L[1,1,t,s] <- 0
      L[1,2,t,s] <- f1.fullarray_MGW$f1[param[s]]
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
    } else {
      L[1,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
      L[1,2,t,s] <- 0
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- 0
    }
  }
}

for(m in 1:M){
  for(s in 1:Sim){
    for(i in 1:I){
      for(k in 1:K){
        N[i,1,k,1,s,m] <- starting_N[i,k]
      }
    }
    for(t in 1:T) {
      for(i in 1:I){
        for(k in 1:K){
          for(j in 1:J){
            if(t %in% removal_index_list[[m]]){
              Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
            } else {
              Y[i,j,k,t,s,m] <- 0
            }
            if(j < J) {
              N[i,j+1,k,t,s,m] <- N[i,j,k,t,s,m] - Y[i,j,k,t,s,m]
            }
          } # end of J for loop
          R[i,k,t,s,m] <- N[i,J,k,t,s,m] - Y[i,J,k,t,s,m] 
        } # end of K for loop
        if(t<T) { # Not including t = T here because with last time step there isn't
          # population growth and movement
          for(k in 1:K){
            # Population growth
            D[i,k,t,s,m] <- trunc(R[i,1,t,s,m]*L[k,1,t,s] + R[i,2,t,s,m]*L[k,2,t,s])
            N[i,1,k,t+1,s,m] <- D[i,k,t,s,m]
          } # end of k for loop
        } # end of if t<T
      } # end of I for loop
    } # end of T for loop
  } # end of Sim for loop
} # end of M for loop


# set up data frame with simulation output
# sum matrices
N_array <- array(NA, dim = c(T,Sim,M))
Y_array <- array(NA, dim = c(T,Sim,M))

for(m in 1:M){
  for(s in 1:Sim){
    for(t in 1:T){
      N_array[t,s,m] <- sum(N[,1,,t,s,m])
      Y_array[t,s,m] <- sum(Y[,1:3,,t,s,m])
    }
  }
}

# Just fish at T=42
# use N_array, not N_summary b/c I need all the simulations for the boxplots
N_matrix_t42 <- N_array[42,,]

M_col2 <- vector()
for(m in 1:M){
  M_col2 <- c(M_col2, rep(m, times = Sim))
}

N_t42 <- vector()
for(m in 1:M){
  for(s in 1:Sim){
    N_t42 <- c(N_t42, N_matrix_t42[s,m])
  }
}

N_t42_df5 <- tibble(M_col2, N_t42)
N_t42_df5$M_col2 <- factor(N_t42_df5$M_col2, levels = seq(1,8,1))


# Saving first two temporal together
# N_t42_df4 and N_t42_df5
N_t42_df4$sim <- "Subset"
N_t42_df5$sim <- "Subset Pre-Spawn"
# sim_output_temporal_fewer <- bind_rows(N_t42_df4, N_t42_df5)


## Temporal, Shorter (Fig 7) ----
I <- 26 # number of pools
J <- 6 # number of passes/gear types
K <- 2 # number of size classes
T <- 42 # removal number (years)
Sim <- 1000 # number of simulations
M <- 9

# random selection of 1000 iterations from 30,000
param <- rdunif(1000,1,30000)

# Starting arrays and values
Y <- array(0, dim = c(I,J,K,T,Sim,M))
N <- array(0, dim = c(I,J,K,T,Sim,M))
R <- array(0, dim = c(I,K,T,Sim,M))
D <- array(0, dim = c(I,K,T,Sim,M))
L <- array(0, dim = c(K,K,T,Sim)) # transition matrix
p <- output_p_MGW # probability of capture dim = c(26,6,42)

# Dates
# original object: dates_MGW
days_btwn <- date_diff_MGW # days between removals
T_hatch <- c(11, 23, 32, 41)

removal_index_vec1 <- vector() # set first type to no date
removal_index_vec2 <- c(1,2) # (2) 3 months
removal_index_vec3 <- seq(1,5) # (5) 6 months
removal_index_vec4 <- seq(1,13) # (13) 12 months
removal_index_vec5 <- seq(1,18) # (18) 18 months
removal_index_vec6 <- seq(1,24) # (24) 24 months
removal_index_vec7 <- seq(1,28) # (28) 30 months
removal_index_vec8 <- seq(1,34) # (34) 36 months
removal_index_vec9 <- seq(1,42) # (41) 47 months

# create list
removal_index_list <- list(removal_index_vec1, removal_index_vec2, removal_index_vec3,
                           removal_index_vec4, removal_index_vec5, removal_index_vec6,
                           removal_index_vec7, removal_index_vec8, removal_index_vec9)

for(s in 1:Sim){
  for(t in 1:T) {
    if(t %in% T_hatch){  # these time points allowing new size 1 individuals 
      L[1,1,t,s] <- 0
      L[1,2,t,s] <- f1.fullarray_MGW$f1[param[s]]
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
    } else {
      L[1,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
      L[1,2,t,s] <- 0
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- 0
    }
  }
}

for(m in 1:M){
  for(s in 1:Sim){
    for(i in 1:I){
      for(k in 1:K){
        N[i,1,k,1,s,m] <- starting_N[i,k]
      }
    }
    for(t in 1:T) {
      for(i in 1:I){
        for(k in 1:K){
          for(j in 1:J){
            if(t %in% removal_index_list[[m]]){
              Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
            } else {
              Y[i,j,k,t,s,m] <- 0
            }
            if(j < J) {
              N[i,j+1,k,t,s,m] <- N[i,j,k,t,s,m] - Y[i,j,k,t,s,m]
            }
          } # end of J for loop
          R[i,k,t,s,m] <- N[i,J,k,t,s,m] - Y[i,J,k,t,s,m] 
        } # end of K for loop
        if(t<T) { # Not including t = T here because with last time step there isn't
          # population growth and movement
          for(k in 1:K){
            # Population growth
            D[i,k,t,s,m] <- trunc(R[i,1,t,s,m]*L[k,1,t,s] + R[i,2,t,s,m]*L[k,2,t,s])
            N[i,1,k,t+1,s,m] <- D[i,k,t,s,m]
          } # end of k for loop
        } # end of if t<T
      } # end of I for loop
    } # end of T for loop
  } # end of Sim for loop
} # end of M for loop

# set up data frame with simulation output
# sum matrices
N_array <- array(NA, dim = c(T,Sim,M))
Y_array <- array(NA, dim = c(T,Sim,M))

for(m in 1:M){
  for(s in 1:Sim){
    for(t in 1:T){
      N_array[t,s,m] <- sum(N[,1,,t,s,m])
      Y_array[t,s,m] <- sum(Y[,1:3,,t,s,m])
    }
  }
}

# Just fish at T=42
# use N_array, not N_summary b/c I need all the simulations for the boxplots
N_matrix_t42 <- N_array[42,,]

M_col2 <- vector()
for(m in 1:M){
  M_col2 <- c(M_col2, rep(m, times = Sim))
}

N_t42 <- vector()
for(m in 1:M){
  for(s in 1:Sim){
    N_t42 <- c(N_t42, N_matrix_t42[s,m])
  }
}


N_t42_df6 <- tibble(M_col2, N_t42) %>%
  filter(M_col2 <=9)
N_t42_df6$M_col2 <- factor(N_t42_df6$M_col2, levels = seq(1,9,1))
# sim_output_temporal_shorter <- N_t42_df6

## Dynamic, Catch (Fig 8) ----
# Hit all pools 1x per year (August), then only revisit 10, 20, ... 100% of pools w/ highest catch during August trip
I <- 26 # number of pools
J <- 6 # number of passes
K <- 2 # number of size classes
T <- 42 # removal number (years)
Sim <- 1000 # number of simulations
M <- 12

# random selection of 1000 iterations from 30,000
param <- rdunif(1000,1,30000)

# Starting arrays and values
Y <- array(0, dim = c(I,J,K,T,Sim,M))
Ysum <- array(NA, dim = c(I,T,Sim,M))
pools <- seq(1, 26)
ysum <- rep(NA, 26)
Ysum_df <- tibble(pools, ysum)
mpercent <- c(0,0,seq(10,100,10))
mpools <- round(mpercent/100*I)
N <- array(0, dim = c(I,J,K,T,Sim,M))
R <- array(0, dim = c(I,K,T,Sim,M))
D <- array(0, dim = c(I,K,T,Sim,M))
L <- array(0, dim = c(K,K,T,Sim)) # transition matrix
p <- output_p_MGW # probability of capture dim = c(26,6,42)

# Dates
# original object: dates_MGW
days_btwn <- date_diff_MGW # days between removals
T_hatch <- c(11, 23, 32, 41)
aug_dates <- c(1,13,24,34) # remove from all pools days

for(s in 1:Sim){
  for(t in 1:T) {
    if(t %in% T_hatch){  # these time points allowing new size 1 individuals 
      L[1,1,t,s] <- 0
      L[1,2,t,s] <- f1.fullarray_MGW$f1[param[s]]
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
    } else {
      L[1,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
      L[1,2,t,s] <- 0
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- 0
    }
  }
}

for(m in 1:M){
  for(s in 1:Sim){
    for(i in 1:I){
      for(k in 1:K){
        N[i,1,k,1,s,m] <- starting_N[i,k]
      }
    }
    for(t in 1:T){
      for(i in 1:I){
        for(k in 1:K){
          for(j in 1:J){
            if(m==1){ # no removals
              Y[i,j,k,t,s,m] <- 0
            }else if(m==2){ # only remove 1x per year 
              if(t %in% aug_dates){
                Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
              }
            }else if(m>2){ # remove all pools 1x per year plus 
              # percentage of pools w/ highest catch
              if(t %in% aug_dates){
                Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
              }else{
                if(i %in% Ysum_df2[1:mpools[m],1]$pools){
                  Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
                }
              }
            }
            if(j < J) {
              N[i,j+1,k,t,s,m] <- N[i,j,k,t,s,m] - Y[i,j,k,t,s,m]
            }
          } # end of J for loop
          R[i,k,t,s,m] <- N[i,J,k,t,s,m] - Y[i,J,k,t,s,m] 
        } # end of K for loop & end of removal section
        if(t<T) { # Not including t = T here because with last time step there isn't
          # population growth and movement
          for(k in 1:K){
            # Population growth
            D[i,k,t,s,m] <- trunc(R[i,1,t,s,m]*L[k,1,t,s] + R[i,2,t,s,m]*L[k,2,t,s])
            N[i,1,k,t+1,s,m] <- D[i,k,t,s,m]
          } # end of k for loop
        } # end of if t<T
        Ysum[i,t,s,m] <- sum(Y[i,,,t,s,m])
        Ysum_df[i,2] <- Ysum[i,t,s,m]
      } # end of I for loop
      Ysum_df2 <- arrange(Ysum_df, desc(ysum)) # this Ysum df changes every t
    } # end of T for loop
  } # end of Sim for loop
} # end of M for loop

# set up data frame with simulation output
# sum matrices
N_array <- array(NA, dim = c(T,Sim,M))
Y_array <- array(NA, dim = c(T,Sim,M))

for(m in 1:M){
  for(s in 1:Sim){
    for(t in 1:T){
      N_array[t,s,m] <- sum(N[,1,,t,s,m])
      Y_array[t,s,m] <- sum(Y[,1:3,,t,s,m])
    }
  }
}

# Just fish at T=42
# use N_array, not N_summary b/c I need all the simulations for the boxplots
N_matrix_t42 <- N_array[42,,]

M_col2 <- vector()
for(m in 1:M){
  M_col2 <- c(M_col2, rep(m, times = Sim))
}

N_t42 <- vector()
for(m in 1:M){
  for(s in 1:Sim){
    N_t42 <- c(N_t42, N_matrix_t42[s,m])
  }
}

N_t42_df7 <- tibble(M_col2, N_t42)
N_t42_df7$M_col2 <- factor(N_t42_df7$M_col2, levels = seq(1,12,1))
# sim_output_dynamic_catch <- N_t42_df7

## Dynamic, Zero (Fig A3) ----
# Once you hit zero M number of times in a pool, stop removing
I <- 26 # number of pools
J <- 6 # number of passes
K <- 2 # number of size classes
T <- 42 # removal number (years)
Sim <- 1000 # number of simulations
M <- 4

# random selection of 1000 iterations from 30,000
param <- rdunif(1000,1,30000)

# Starting arrays and values
Y <- array(0, dim = c(I,J,K,T,Sim,M))
Ysum <- array(NA, dim = c(I,T,Sim,M))
N <- array(0, dim = c(I,J,K,T,Sim,M))
R <- array(0, dim = c(I,K,T,Sim,M))
D <- array(0, dim = c(I,K,T,Sim,M))
L <- array(0, dim = c(K,K,T,Sim)) # transition matrix

p <- output_p_MGW # probability of capture dim = c(26,6,42)

# Dates
# original object: dates_MGW
days_btwn <- date_diff_MGW # days between removals
T_hatch <- c(11, 23, 32, 41)

for(s in 1:Sim){
  for(t in 1:T) {
    if(t %in% T_hatch){  # these time points allowing new size 1 individuals 
      L[1,1,t,s] <- 0
      L[1,2,t,s] <- f1.fullarray_MGW$f1[param[s]]
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
    } else {
      L[1,1,t,s] <- s1.fullarray_MGW$s1[param[s]]^days_btwn[t]
      L[1,2,t,s] <- 0
      L[2,2,t,s] <- s2.fullarray_MGW$s2[param[s]]^days_btwn[t]
      L[2,1,t,s] <- 0
    }
  }
}

for(m in 1:M){
  for(s in 1:Sim){
    for(i in 1:I){
      for(k in 1:K){
        N[i,1,k,1,s,m] <- starting_N[i,k]
      }
    }
    for(t in 1:T){
      for(i in 1:I) {
        for(k in 1:K){
          for(j in 1:J) {
            if(m==4){ # keep removing as normal
              Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
            }else if(m==1){ # stop removals after 1 date with zero caught
              if(t > 1){
                if(Ysum[i,t-1,s,m] == 0){
                  Y[i,j,k,t,s,m] <- 0
                }else{
                  Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
                }
              }else{
                Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
              }
            }else if(m==2){ # stop removals after 2 dates with zero caught
              if(t > 2){
                if(Ysum[i,t-1,s,m] == 0 & Ysum[i,t-2,s,m] == 0){
                  Y[i,j,k,t,s,m] <- 0
                }else{
                  Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
                }
              }else{
                Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
              }
            }else if(m==3){ # stop removals after 3 dates with zero caught
              if(t > 3){
                if(Ysum[i,t-1,s,m] == 0 & Ysum[i,t-2,s,m] == 0 & Ysum[i,t-3,s,m] == 0){
                  Y[i,j,k,t,s,m] <- 0
                }else{
                  Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
                }
              }else{
                Y[i,j,k,t,s,m] <- rbinom(1,N[i,j,k,t,s,m],p[i,j,t]) # removal
              }
            }
            if(j < J) {
              N[i,j+1,k,t,s,m] <- N[i,j,k,t,s,m] - Y[i,j,k,t,s,m]
            }
          } # end of J for loop
          R[i,k,t,s,m] <- N[i,J,k,t,s,m] - Y[i,J,k,t,s,m] 
        } # end of K for loop 
        if(t<T) { # Not including t = T here because with last time step there isn't
          # population growth and movement
          for(k in 1:K){
            # Population growth
            D[i,k,t,s,m] <- trunc(R[i,1,t,s,m]*L[k,1,t,s] + R[i,2,t,s,m]*L[k,2,t,s])
            N[i,1,k,t+1,s,m] <- D[i,k,t,s,m]
          } # end of k for loop
        } # end of if t<T
        Ysum[i,t,s,m] <- sum(Y[i,,,t,s,m])
      } # end of I for loop
    } # end of T for loop
  } # end of Sim for loop
} # end of M for loop

# set up data frame with simulation output
# sum matrices
N_array <- array(NA, dim = c(T,Sim,M))
Y_array <- array(NA, dim = c(T,Sim,M))

for(m in 1:M){
  for(s in 1:Sim){
    for(t in 1:T){
      N_array[t,s,m] <- sum(N[,1,,t,s,m])
      Y_array[t,s,m] <- sum(Y[,1:3,,t,s,m])
    }
  }
}

# Just fish at T=42
# use N_array, not N_summary b/c I need all the simulations for the boxplots
N_matrix_t42 <- N_array[42,,]

M_col2 <- vector()
for(m in 1:M){
  M_col2 <- c(M_col2, rep(m, times = Sim))
}

N_t42 <- vector()
for(m in 1:M){
  for(s in 1:Sim){
    N_t42 <- c(N_t42, N_matrix_t42[s,m])
  }
}

N_t42_df8 <- tibble(M_col2, N_t42)
N_t42_df8$M_col2 <- factor(N_t42_df8$M_col2, levels = seq(1,4,1))
# sim_output_dynamic_zero <- N_t42_df8

# Plotting ----
## Figure 2 Bar plot ----
# EAC original code
# bar plot
p1 <- ggplot(data = barplot_data_eac, aes(x = date, y = count, fill = sizeclass)) + 
  geom_bar(stat = "identity", position = "stack", width = 8) +
  theme_classic() +
  labs(y = "", x = "Time") +
  scale_fill_discrete(name = "Size Class (mm)", type = c("#1f78b4","#a6cee3")) +
  theme(legend.position = "none", text = element_text(size=12)) +
  scale_x_continuous(breaks = as.Date(c("2014-01-01", "2015-01-01")), 
                     labels = c("Jan 2014", "Jan 2015")) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 399)) + 
  geom_segment(aes(x = as.Date("2014-08-18"), y = 30, 
                   xend = as.Date("2014-08-18"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2014-11-26"), y = 30, 
                   xend = as.Date("2014-11-26"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2015-01-21"), y = 30, 
                   xend = as.Date("2015-01-21"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm")))+ 
  geom_segment(aes(x = as.Date("2015-02-25"), y = 30, 
                   xend = as.Date("2015-02-25"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2015-04-13"), y = 30, 
                   xend = as.Date("2015-04-13"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2015-05-27"), y = 30, 
                   xend = as.Date("2015-05-27"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm")))


# MGW original code
# bar plot
p2 <- ggplot(data = barplot_data_mgw, aes(x = date, y = count, fill = sizeclass)) + 
  theme_classic() +
  geom_bar(stat = "identity", position = "stack", width = 8) +
  labs(y = "Green Sunfish (#)", x = "Time") +
  scale_fill_discrete(name = "Size Class (mm)", type = c("#1f78b4","#a6cee3")) +
  theme(legend.position = c(0.8, 0.8), text = element_text(size=12),
                      labels = c("< 50", "\u2265 50")) +
  scale_x_continuous(breaks = as.Date(c("2018-01-01", "2019-01-01", 
                                        "2020-01-01", "2021-01-01")), 
                     labels = c("Jan 2018", "Jan 2019", 
                                "Jan 2020", "Jan 2021")) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 2200)) + 
  geom_segment(aes(x = as.Date("2019-07-22"), y = 150, 
                   xend = as.Date("2019-07-22"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2020-02-27"), y = 150, 
                   xend = as.Date("2020-02-27"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2020-03-27"), y = 150, 
                   xend = as.Date("2020-03-27"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2020-05-20"), y = 150, 
                   xend = as.Date("2020-05-20"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2020-11-13"), y = 150, 
                   xend = as.Date("2020-11-13"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2020-12-16"), y = 150, 
                   xend = as.Date("2020-12-16"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2021-04-01"), y = 150, 
                   xend = as.Date("2021-04-01"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2021-04-28"), y = 150, 
                   xend = as.Date("2021-04-28"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = as.Date("2021-07-07"), y = 150, 
                   xend = as.Date("2021-07-07"), yend = 0),
               arrow = arrow(length = unit(0.15, "cm"))) 

plot_grid(p2, p1, labels = c("A", "B"), nrow =1, rel_widths = c(50,23))

ggsave("Figure-2.png", units = "in", width = 8, height = 3, dpi = 300)


## Figure 4 Model results ----
# McGee Wash
N.sum1 <- tibble(N.sum1 = N.sum.fullarray_MGW$N.sum[,1])
p1 <- ggplot(data = N.sum1) + 
  geom_density(aes(x = N.sum1), fill = "#08519c", alpha = 0.5, color = "#08519c") + 
  theme_classic() + 
  xlab("Green Sunfish (#)") + ylab("Density") + 
  theme(text = element_text(size = 7))
N.sum42 <- tibble(N.sum42 = N.sum.fullarray_MGW$N.sum[,42])
p2 <- ggplot(data = N.sum42) + 
  geom_density(aes(x = N.sum42), fill = "#08519c", alpha = 0.5, color = "#08519c") + 
  theme_classic() + 
  xlab("Green Sunfish (#)") + ylab("Density") + 
  theme(text = element_text(size = 7))
# Main Plot
p3 <- ggplot(pop_time_MGW, aes(x = dates, y = Y)) +
  geom_bar(stat = "identity", fill = "#525252", width = 15) + 
  xlab("Time") + ylab("Green Sunfish (#)") +
  geom_ribbon(aes(x = dates, y = N, ymin = CImin, ymax = CImax), 
              fill = "grey70") +
  geom_line(aes(x = dates, y = N), color = "#08519c", size = 0.5) +
  #annotate(geom = "text", x = as_date("2018-03-25"), y = 5500, label = "Total population", color = "#08519c", size = 3) +
  #annotate(geom = "text", x = as_date("2017-11-20"), y = 2500, label = "Fish\nremoved", color = "black", size = 3) + 
  theme_classic() + 
  theme(text = element_text(size = 10)) + 
  scale_x_continuous(breaks = as.Date(c("2017-08-09", "2018-01-01", "2019-01-01", 
                                        "2020-01-01", "2021-01-01", "2021-07-07")),
                     labels = c("Aug 8\n2017", "Jan 1\n2018", "Jan 1\n2019", 
                                "Jan 1\n2020", "Jan 1\n2021", "Jul 7\n2021")) + 
  scale_y_continuous(expand = c(0,0)) 

# Combine Plots
p4 <- ggdraw(p3) + draw_plot(p1, x = 0.19, y = 0.69, width = 0.32, height = 0.32) + 
  draw_plot(p2,  x = 0.66, y = 0.22, width = 0.32, height = 0.32)


# East Ash Creek
N.sum1 <- tibble(N.sum1 = N.sum.fullarray_EAC$N.sum[,1])
p1 <- ggplot(data = N.sum1) + 
  geom_density(aes(x = N.sum1), fill = "#08519c", alpha = 0.5, color = "#08519c") + 
  theme_classic() + 
  xlab("Green Sunfish (#)") + ylab("Density") + 
  theme(text = element_text(size = 7))
N.sum17 <- tibble(N.sum17 = N.sum.fullarray_EAC$N.sum[,17])
p2 <- ggplot(data = N.sum17) + 
  geom_density(aes(x = N.sum17), fill = "#08519c", alpha = 0.5, color = "#08519c") + 
  theme_classic() + 
  xlab("Green Sunfish (#)") + ylab("Density") + 
  theme(text = element_text(size = 7))
p3 <- ggplot(pop_time_EAC, aes(x = dates, y = Y)) +
  geom_bar(stat = "identity", fill = "#525252", width = 15) + 
  xlab("Time") + ylab("Green Sunfish (#)") +
  geom_ribbon(aes(x = dates, y = N, ymin = CImin, ymax = CImax), 
              fill = "grey70") +
  geom_line(aes(x = dates, y = N), color = "#08519c", size = 0.5) +
  #annotate(geom = "text", x = as_date("2018-03-25"), y = 5500, label = "Total population", color = "#08519c", size = 3) +
  #annotate(geom = "text", x = as_date("2017-11-20"), y = 2500, label = "Fish\nremoved", color = "black", size = 3) + 
  theme_classic() + 
  theme(text = element_text(size = 10)) + 
  scale_x_continuous(breaks = as.Date(c("2013-11-05", "2014-01-01", "2014-07-01", "2015-01-01", "2015-05-27")),
                     labels = c("Nov 5\n2013", "Jan 1\n2014","Jul 1\n2014", "Jan 1\n2015", "May 27\n2015")) + 
  scale_y_continuous(expand = c(0,0)) 

p5 <- ggdraw(p3) + draw_plot(p1, x = 0.275, y = 0.69, width = 0.32, height = 0.32) + 
  draw_plot(p2,  x = 0.66, y = 0.20, width = 0.32, height = 0.32)


# plotting MGW and EAC together in a multipanel plot
plot_grid(p4, p5, labels = c("A", "B"), nrow =2)

ggsave("Figure-4.png", units = "in", width = 5, height = 6, dpi = 300)


## Simulation plots ----
### Figure 5 ----
# subset w/out Mcol = 1
sim_output_spatial_random <- sim_output_spatial_random %>% filter(M_col2 != 1)
sim_output_spatial_random$M_col2 <- factor(sim_output_spatial_random$M_col2, levels = seq(2,11,1))

# subset of Random
random_subset <- sim_output_spatial_random %>% filter(sim == "Random")
M_col2 <- seq(2,11,1)
sim <- rep("Random", 10)
q2.5 <- rep(NA, 10)
q25 <- rep(NA, 10)
q50 <- rep(NA, 10)
q75 <- rep(NA, 10)
q97.5 <- rep(NA,10)
random_summary <- tibble(M_col2, sim, q2.5, q25, q50, q75, q97.5)
random_summary$M_col2 <- factor(random_summary$M_col2, levels = seq(2,11,1))

for(i in 2:11){
  subset <- random_subset %>% filter(M_col2 == i)
  random_summary[random_summary$M_col2 == i,]$q2.5 <- quantile(subset$N_t42, prob = 0.025)
  random_summary[random_summary$M_col2 == i,]$q25 <- quantile(subset$N_t42, prob = 0.25)
  random_summary[random_summary$M_col2 == i,]$q50 <- quantile(subset$N_t42, prob = 0.50)
  random_summary[random_summary$M_col2 == i,]$q75 <- quantile(subset$N_t42, prob = 0.75)
  random_summary[random_summary$M_col2 == i,]$q97.5 <- quantile(subset$N_t42, prob = 0.975)
}

# subset of Contiguous
contig_subset <- sim_output_spatial_random %>% filter(sim == "Random Contiguous")
M_col2 <- seq(2,11,1)
sim <- rep("Random Contiguous", 10)
q2.5 <- rep(NA, 10)
q25 <- rep(NA, 10)
q50 <- rep(NA, 10)
q75 <- rep(NA, 10)
q97.5 <- rep(NA,10)
contig_summary <- tibble(M_col2, sim, q2.5, q25, q50, q75, q97.5)
contig_summary$M_col2 <- factor(contig_summary$M_col2, levels = seq(2,11,1))

for(i in 2:11){
  subset <- contig_subset %>% filter(M_col2 == i)
  contig_summary[contig_summary$M_col2 == i,]$q2.5 <- quantile(subset$N_t42, prob = 0.025)
  contig_summary[contig_summary$M_col2 == i,]$q25 <- quantile(subset$N_t42, prob = 0.25)
  contig_summary[contig_summary$M_col2 == i,]$q50 <- quantile(subset$N_t42, prob = 0.50)
  contig_summary[contig_summary$M_col2 == i,]$q75 <- quantile(subset$N_t42, prob = 0.75)
  contig_summary[contig_summary$M_col2 == i,]$q97.5 <- quantile(subset$N_t42, prob = 0.975)
}

ggplot(data = random_summary, aes(x = M_col2, y = q50)) +
  geom_boxplot(stat = "identity",
               aes(lower = q25, upper = q75, middle = q50, 
                   ymin = q2.5, ymax = q97.5, group = M_col2), 
               color = "grey15", fill = "#8856a7", width = 0.7) +
  geom_boxplot(data = contig_summary, stat = "identity",
               aes(lower = q25, upper = q75, middle = q50, 
                   ymin = q2.5, ymax = q97.5, group = M_col2), 
               color = "#1b9e77", alpha = 0.4, width = 0.7) +
  theme_bw() + 
  ylab("Final Green Sunfish Population (#)") + 
  xlab("Simulated Pools Chosen for Removal (%)") + 
  scale_x_discrete(breaks = seq(2,11,1), labels = seq(10,100,10)) +
  geom_hline(yintercept = sum(starting_N), lty = 2) +
  theme(text = element_text(size = 12)) + 
  annotate(geom = "text", x = 1.5, y = 70, label = "Random Contiguous Pool Selection", size = 4, hjust = 0, color = "#1b9e77") + 
  annotate(geom = "text", x = 1.5, y = 140, label = "Random Pool Selection", size = 4, hjust = 0, color = "#8856a7") + 
  annotate(geom = "text", x = 1.1, y = 11000, label = "Starting\nPopulation", size = 3, hjust = 0)+
  scale_y_continuous(trans = scales::pseudo_log_trans(),labels = label_comma(),
                     breaks = c(0, 10,100, 1000, 10000, 100000))

ggsave("Figure-5.png", units = "in", width = 6, height = 4, dpi = 300)

### Figure 6 ----
sim_output_temporal_fewer$M_col2 <- factor(sim_output_temporal_fewer$M_col2, levels = seq(1,8,1))

# plain subset
plain_subset <- sim_output_temporal_fewer %>% filter(sim == "Subset")
M_col2 <- seq(1,8,1)
sim <- rep("Subset", 8)
q2.5 <- rep(NA, 8)
q25 <- rep(NA, 8)
q50 <- rep(NA, 8)
q75 <- rep(NA, 8)
q97.5 <- rep(NA,8)
plain_summary <- tibble(M_col2, sim, q2.5, q25, q50, q75, q97.5)
plain_summary$M_col2 <- factor(plain_summary$M_col2, levels = seq(1,8,1))

for(i in 1:8){
  subset <- plain_subset %>% filter(M_col2 == i)
  plain_summary[plain_summary$M_col2 == i,]$q2.5 <- quantile(subset$N_t42, prob = 0.025)
  plain_summary[plain_summary$M_col2 == i,]$q25 <- quantile(subset$N_t42, prob = 0.25)
  plain_summary[plain_summary$M_col2 == i,]$q50 <- quantile(subset$N_t42, prob = 0.50)
  plain_summary[plain_summary$M_col2 == i,]$q75 <- quantile(subset$N_t42, prob = 0.75)
  plain_summary[plain_summary$M_col2 == i,]$q97.5 <- quantile(subset$N_t42, prob = 0.975)
}

# pre-spawn subset
prespawn_subset <- sim_output_temporal_fewer %>% filter(sim == "Subset Pre-Spawn")
M_col2 <- seq(1,8,1)
sim <- rep("Subset Pre-Spawn", 8)
q2.5 <- rep(NA, 8)
q25 <- rep(NA, 8)
q50 <- rep(NA, 8)
q75 <- rep(NA, 8)
q97.5 <- rep(NA,8)
prespawn_summary <- tibble(M_col2, sim, q2.5, q25, q50, q75, q97.5)
prespawn_summary$M_col2 <- factor(prespawn_summary$M_col2, levels = seq(1,8,1))

for(i in 1:8){
  subset <- prespawn_subset %>% filter(M_col2 == i)
  prespawn_summary[prespawn_summary$M_col2 == i,]$q2.5 <- quantile(subset$N_t42, prob = 0.025)
  prespawn_summary[prespawn_summary$M_col2 == i,]$q25 <- quantile(subset$N_t42, prob = 0.25)
  prespawn_summary[prespawn_summary$M_col2 == i,]$q50 <- quantile(subset$N_t42, prob = 0.50)
  prespawn_summary[prespawn_summary$M_col2 == i,]$q75 <- quantile(subset$N_t42, prob = 0.75)
  prespawn_summary[prespawn_summary$M_col2 == i,]$q97.5 <- quantile(subset$N_t42, prob = 0.975)
}

ggplot(data = plain_summary, aes(x = M_col2, y = q50)) +
  geom_boxplot(stat = "identity",
               aes(lower = q25, upper = q75, middle = q50, 
                   ymin = q2.5, ymax = q97.5, group = M_col2), 
               color = "grey15", fill = "#8856a7") + 
  geom_boxplot(data = prespawn_summary, stat = "identity", 
               aes(lower = q25, upper = q75, middle = q50, 
                   ymin = q2.5, ymax = q97.5, group = M_col2), 
               color = "#1b9e77", alpha = 0.4) +
  theme_bw() + 
  ylab("Final Green Sunfish Population (#)") + 
  xlab("Simulated Removal Events (#)") + 
  scale_x_discrete(breaks = seq(1,8,1), labels = c(5,10,15,20,25,30,35,42)) +
  geom_hline(yintercept = sum(starting_N), lty = 2) +
  theme(text = element_text(size = 12)) + 
  annotate(geom = "text", x = 4, y = 150000, label = "Subset of Removal Events", size = 4, hjust = 0, color = "#8856a7") + 
  annotate(geom = "text", x = 4, y = 75000, label = "Subset Prioritizing Pre-spawn Removal", size = 4, hjust = 0, color = "#1b9e77") + 
  annotate(geom = "text", x = 1, y = 11000, label = "Starting\nPopulation", size = 3, hjust = 0)+
  scale_y_continuous(trans = scales::pseudo_log_trans(),labels = label_comma(),
                     breaks = c(0, 10,100, 1000, 10000, 100000))

ggsave("Figure-6.png", units = "in", width = 6, height = 4, dpi = 300)


### Figure 7 ----
# subset without M_col2 = 1
sim_output_temporal_shorter <- sim_output_temporal_shorter %>% filter(M_col2 != 1)
sim_output_temporal_shorter$M_col2 <- factor(sim_output_temporal_shorter$M_col2, levels = seq(2,9,1))


M_col2 <- seq(2,9,1)
q2.5 <- rep(NA, 8)
q25 <- rep(NA, 8)
q50 <- rep(NA, 8)
q75 <- rep(NA, 8)
q97.5 <- rep(NA,8)
summary <- tibble(M_col2, q2.5, q25, q50, q75, q97.5)
summary$M_col2 <- factor(summary$M_col2, levels = seq(2,9,1))

for(i in 2:9){
  subset <- sim_output_temporal_shorter %>% filter(M_col2 == i)
  summary[summary$M_col2 == i,]$q2.5 <- quantile(subset$N_t42, prob = 0.025)
  summary[summary$M_col2 == i,]$q25 <- quantile(subset$N_t42, prob = 0.25)
  summary[summary$M_col2 == i,]$q50 <- quantile(subset$N_t42, prob = 0.50)
  summary[summary$M_col2 == i,]$q75 <- quantile(subset$N_t42, prob = 0.75)
  summary[summary$M_col2 == i,]$q97.5 <- quantile(subset$N_t42, prob = 0.975)
}


ggplot(data = summary, aes(x = M_col2, y = q50)) +
  geom_boxplot(stat = "identity",
               aes(lower = q25, upper = q75, middle = q50, 
                   ymin = q2.5, ymax = q97.5, group = M_col2)) + 
  theme_bw() + 
  ylab("Final Green Sunfish Population (#)") + 
  xlab("Duration of Simulated Removal Program (# months)") + 
  geom_hline(yintercept = sum(starting_N), lty = 2) + 
  scale_x_discrete(breaks = seq(2,9,1), 
                   labels = c(3,6,12,18,24,30,36,47)) +
  theme(text = element_text(size = 12)) +
  annotate(geom = "text", x = 1, y = 11000, label = "Starting\nPopulation", size = 3, hjust = 0) +
  scale_y_continuous(trans = scales::pseudo_log_trans(),labels = label_comma(),
                     breaks = c(0, 10,100, 1000, 10000, 100000))

ggsave("Figure-7.png", units = "in", width = 6, height = 4, dpi = 300)


### Figure 8 ----
# subset without M_col2 = 1
sim_output_dynamic_catch <- sim_output_dynamic_catch %>% filter(M_col2 != 1)
sim_output_dynamic_catch$M_col2 <- factor(sim_output_dynamic_catch$M_col2, levels = seq(2,12,1))


M_col2 <- seq(2,12,1)
q2.5 <- rep(NA, 11)
q25 <- rep(NA, 11)
q50 <- rep(NA, 11)
q75 <- rep(NA, 11)
q97.5 <- rep(NA,11)
summary <- tibble(M_col2, q2.5, q25, q50, q75, q97.5)
summary$M_col2 <- factor(summary$M_col2, levels = seq(2,12,1))

for(i in 2:12){
  subset <- sim_output_dynamic_catch %>% filter(M_col2 == i)
  summary[summary$M_col2 == i,]$q2.5 <- quantile(subset$N_t42, prob = 0.025)
  summary[summary$M_col2 == i,]$q25 <- quantile(subset$N_t42, prob = 0.25)
  summary[summary$M_col2 == i,]$q50 <- quantile(subset$N_t42, prob = 0.50)
  summary[summary$M_col2 == i,]$q75 <- quantile(subset$N_t42, prob = 0.75)
  summary[summary$M_col2 == i,]$q97.5 <- quantile(subset$N_t42, prob = 0.975)
}

ggplot(data = summary, aes(x = M_col2, y = q50)) +
  geom_boxplot(stat = "identity",
               aes(lower = q25, upper = q75, middle = q50, 
                   ymin = q2.5, ymax = q97.5, group = M_col2)) + 
  theme_bw() + 
  ylab("Final Green Sunfish Population (#)") + 
  xlab("Removals from Pools with Highest Catch (%)") + 
  geom_hline(yintercept = sum(starting_N), lty = 2) + 
  scale_x_discrete(breaks = seq(2,12,1), 
                   labels = c("All pools\n once/year only", "10", "20", "30",
                              "40", "50", "60", "70", "80", "90", "All pools")) +
  theme(text = element_text(size = 12)) +
  annotate(geom = "text", x = 0.7, y = 11000, label = "Starting\nPopulation", size = 3, hjust = 0) +
  scale_y_continuous(trans = scales::pseudo_log_trans(),labels = label_comma(),
                     breaks = c(0, 10,100, 1000, 10000, 100000))

ggsave("Figure-8.png", units = "in", width = 7, height = 5, dpi = 300)


### Figure A2 ----
# subset without Mcol = 1
sim_output_spatial_vol <- sim_output_spatial_vol %>% filter(M_col2 != 1)
sim_output_spatial_vol$M_col2 <- factor(sim_output_spatial_vol$M_col2, levels = seq(2,11,1))

M_col2 <- seq(2,11,1)
q2.5 <- rep(NA, 10)
q25 <- rep(NA, 10)
q50 <- rep(NA, 10)
q75 <- rep(NA, 10)
q97.5 <- rep(NA,10)
summary <- tibble(M_col2, q2.5, q25, q50, q75, q97.5)
summary$M_col2 <- factor(summary$M_col2, levels = seq(2,11,1))

for(i in 2:11){
  subset <- sim_output_spatial_vol %>% filter(M_col2 == i)
  summary[summary$M_col2 == i,]$q2.5 <- quantile(subset$N_t42, prob = 0.025)
  summary[summary$M_col2 == i,]$q25 <- quantile(subset$N_t42, prob = 0.25)
  summary[summary$M_col2 == i,]$q50 <- quantile(subset$N_t42, prob = 0.50)
  summary[summary$M_col2 == i,]$q75 <- quantile(subset$N_t42, prob = 0.75)
  summary[summary$M_col2 == i,]$q97.5 <- quantile(subset$N_t42, prob = 0.975)
}

ggplot(data = summary, aes(x = M_col2, y = q50)) +
  geom_boxplot(stat = "identity",
               aes(lower = q25, upper = q75, middle = q50, 
                   ymin = q2.5, ymax = q97.5, group = M_col2)) + 
  theme_minimal() + 
  ylab("Final Green Sunfish Population (#)") + 
  xlab("Simulated Pools Chosen for Removal (Largest % by Volume)") + 
  geom_hline(yintercept = sum(starting_N), lty = 2) + 
  scale_x_discrete(breaks = seq(2,11,1), labels = seq(10,100,10)) + 
  annotate(geom = "text", x = 1, y = 35000, label = "Starting\nPopulation", size = 3, hjust = 0) +
  scale_y_continuous(labels = label_comma())

ggsave("Figure-A2.png", units = "in", width = 6, height = 4, dpi = 300)

### Figure A3 ----
sim_output_dynamic_zero$M_col2 <- factor(sim_output_dynamic_zero$M_col2, levels = seq(1,4,1))

M_col2 <- seq(1,4,1)
q2.5 <- rep(NA, 4)
q25 <- rep(NA, 4)
q50 <- rep(NA, 4)
q75 <- rep(NA, 4)
q97.5 <- rep(NA,4)
summary <- tibble(M_col2, q2.5, q25, q50, q75, q97.5)
summary$M_col2 <- factor(summary$M_col2, levels = seq(1,4,1))

for(i in 1:4){
  subset <- sim_output_dynamic_zero %>% filter(M_col2 == i)
  summary[summary$M_col2 == i,]$q2.5 <- quantile(subset$N_t42, prob = 0.025)
  summary[summary$M_col2 == i,]$q25 <- quantile(subset$N_t42, prob = 0.25)
  summary[summary$M_col2 == i,]$q50 <- quantile(subset$N_t42, prob = 0.50)
  summary[summary$M_col2 == i,]$q75 <- quantile(subset$N_t42, prob = 0.75)
  summary[summary$M_col2 == i,]$q97.5 <- quantile(subset$N_t42, prob = 0.975)
}

ggplot(data = summary, aes(x = M_col2, y = q50)) +
  geom_boxplot(stat = "identity",
               aes(lower = q25, upper = q75, middle = q50, 
                   ymin = q2.5, ymax = q97.5, group = M_col2)) + 
  theme_minimal() + 
  ylab("Final Green Sunfish Population (#)") + 
  xlab("Ceasing Sampling in Pools") + 
  scale_x_discrete(breaks = seq(1,4,1), 
                   labels = c("Zero Fish 1x", 
                              "Zero Fish 2x", 
                              "Zero Fish 3x", 
                              "Continued removals")) +
  theme(text = element_text(size = 12)) 

ggsave("Figure-A3.png", units = "in", width = 7, height = 5, dpi = 300)

