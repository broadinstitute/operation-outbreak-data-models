args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  data_folder <- "."
} else {
  data_folder <- args[1]  
}

library(dplyr)
set.seed(594709947L)
library(ggplot2)
theme_set(theme_bw())
library(plyr)
library(reshape2)
library(magrittr)
library(pomp)
stopifnot(packageVersion("pomp")>="2.0")

# INPUT ================================================================

input_params <- read.csv(file.path(data_folder, "input_params.csv"))

# Duration of the outbreak in minutes
total_min <- input_params$T
total_qrt <- total_min / 15

# Total number of participants, and number of initial cases
pop_size <- input_params$N
init_cases <- input_params$I0

# Expected cases
case_data <- read.csv(file.path(data_folder, "cases.csv"))
case_data <- data.frame(time = case_data$Time, cases = case_data$Count)
#head(case_data)

# CLEAN DATA ================================================================

# convert minute data to quarter data 
case_data$quarter = floor(case_data$time / 15)

qrt_data <- data.frame(quarter = case_data$quarter, cases = case_data$cases)
qrt_data <- group_by(qrt_data, quarter)

sum_cases = c()
for (i in unique(qrt_data$quarter)) {
  temp = subset(qrt_data, quarter == i)
  sum_cases = c(sum_cases, sum(temp$cases))
}

qrt_data <- data.frame(quarter = unique(qrt_data$quarter), cases = sum_cases)

other_data = data.frame('time' = seq(1, max(total_min, pop_size)), 'cases' = rep(0, max(total_min, pop_size)))

# keep the relevant columns 
case_data = case_data[,1:2]

#FINAL data
temp = merge(x = other_data, y = case_data, by = "time", all.x = TRUE)
temp[is.na(temp)] <- 0
data = data.frame('time' = temp$time, 'cases' = temp$cases.x + temp$cases.y)
#data

# plot both by minute and quarter-of-hour curves 
ggplot(data=subset(data, time <= total_min), aes(x=time, y=cases, group=1)) + geom_line()
ggsave(file.path(data_folder, "minute.pdf"))

ggplot(data=subset(qrt_data, quarter <= total_qrt), aes(x=quarter, y=cases, group=1)) + geom_line()
ggsave(file.path(data_folder, "quarter.pdf"))

data=subset(data, time <= total_min)

# COMPARTMENTAL MODEL ===========================================================

# The C snippets that define the dynamics of the system
rproc <- Csnippet("
                  double beta, foi;
                  double rate[2], trans[2];
                  
                  // transmission rate
                  beta = R0 * gamma;
                  
                  // expected force of infection
                  foi = beta * I/pop;
                  
                  rate[0] = foi;       // stochastic force of infection
                  rate[1] = gamma;     // recovery rate
                  
                  // transitions between classes
                  reulermultinom(1, S, &rate[0], dt, &trans[0]);
                  reulermultinom(1, I, &rate[1], dt, &trans[1]);
                  
                  C = trans[0];            // true number of new infectious cases in a single day
                  D = trans[1];            // true number of cases that are removed from the population in a single day
                  
                  S = S - C;
                  I = I + C - D;
                  
                  // Population is constant, so we can calculate R like so (instead of R = R + D):
                  R = pop - S - I;
                  ")

initlz <- Csnippet("
                   double m = pop/(S_0 + I_0 + R_0);
                   
                   S = nearbyint(m*S_0);
                   I = nearbyint(m*I_0);
                   R = nearbyint(m*R_0);
                   
                   C = 0;
                   D = 0;
                   ")

dmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  ")

rmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")

i0 <- init_cases / pop_size
s0 <- 1 - i0
fixed_params <- c(pop = pop_size, S_0 = s0, I_0 = i0, R_0 = 0)
fixed_names <- c("pop", "S_0", "I_0", "R_0")

# random parameters that will estimate
rp_names <- c("R0", "gamma", "rho", "psi")

all_param_names=c(rp_names, fixed_names)

data %>% 
  pomp(t0 = data$time[1],
       time = "time",
       rprocess = euler(rproc, delta.t=1),
       rinit = initlz,
       dmeasure = dmeas,
       rmeasure = rmeas,
       accumvars = c("C", "D"),
       statenames = c("S", "I", "R", "C", "D"),
       partrans = parameter_trans(
         log = c("psi","R0"),
         logit = c("rho", "gamma"),
         barycentric = c("S_0","I_0","R_0")),
       paramnames = c(rp_names, fixed_names)
  ) -> m1

m1 %>% as.data.frame() %>%
  melt(id="time") %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  facet_grid(variable~., scales = "free_y")

# GLOBAL MLE ================================================================

library(foreach)
library(doParallel)

registerDoParallel()
set.seed(998468235L, kind="L'Ecuyer")  

param_box <- rbind(
  R0 = c(0.1,5),
  gamma = c(0,1),
  rho = c(0.5,1),
  psi = c(0,1)
)

stew(file=file.path(data_folder, "box_search_global_with_intervention.rda"),{
  w2 <- getDoParWorkers()
  t2 <- system.time({
    m2 <- foreach(i=1:20,.packages='pomp',.combine=rbind,
                  .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      guess <- apply(param_box,1,function(x)runif(1,x[1],x[2]))
      mf <- mif2(
        m1,
        params=c(guess, fixed_params),
        Np=5000,
        Nmif=100,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        rw.sd=rw.sd(R0=0.02,gamma=0.02,rho=0.02,psi=0.02)
      )
      ll <- logmeanexp(replicate(5,logLik(pfilter(mf,Np=5000))), se = TRUE)
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
    }
  })
}, seed=290860873, kind="L'Ecuyer")

library(plyr)
global_params <- arrange(rbind(m2),-loglik)
write.csv(global_params, file=file.path(data_folder, "global_params.csv"), row.names=FALSE, na="")

# SIMULATION W/ BEST PARAMETERS =================================================================

log_idx <- length(global_params) - 1
mle_global <- global_params[which.max( global_params[,log_idx] ), ] 
mle_global %>% extract(all_param_names) %>% unlist() -> theta_global

#print(theta_global)
write.csv(theta_global,file=file.path(data_folder, "best_params.csv"),row.names=TRUE,na="")

m1  %>%
  simulate(params=c(theta_global, fixed_params),nsim=9,
           format = "data.frame",include.data=TRUE) -> sim_data

sim_data %>%
  ggplot(aes(x=time,y=cases,group=.id,color=(.id=="data")))+
  guides(color=FALSE)+
  geom_line()+facet_wrap(~.id,ncol=2)

ggsave(file.path(data_folder, "simulations.pdf"))

# COMPARE CUMULATIVE NUMBERS OF OUTBREAK WITH ACTUAL DATA ======================================

num_sim <- 100

m1 %>% 
  simulate(params=c(theta_global), 
           nsim=num_sim,format = "data.frame",include.data=TRUE) -> simulation_data

cumulative_data = c()

# Getting the cumulative curve for observed data
sum_data = 0
for (i in 1:total_min) {
  sum_data = sum_data + data$cases[i]
  cumulative_data = c(cumulative_data, sum_data)
}

# Getting the average cumulative curve for the simulations
total_size = c()
for (i in 1:total_min) {
  temp = subset(simulation_data, .id == i)
  s = sum(temp$C)
  total_size = c(total_size, s)
}

# Taking the median
sel_sim = 0.5 * num_sim
median_index = order(total_size)[sel_sim]

sum_int = 0
cumulative_sim = c()
for (i in 1:total_min) {
  sim_int <- subset(simulation_data, .id == median_index)
  sum_int = sum_int + sim_int$C[i]
  cumulative_sim = c(cumulative_sim, sum_int)
}

cumulative_numbers = data.frame('time' = seq(1,total_min), 
                                'actual_data' = cumulative_data,
                                'simulated_data' = cumulative_sim)

ggplot(cumulative_numbers, aes(time)) + 
  geom_line(aes(y = actual_data, colour = "Real Data")) + 
  geom_line(aes(y = simulated_data, colour = "Simulated")) +
  ylab('Cumulative Number of Cases')

ggsave(file.path(data_folder, "cumulative.pdf"))
