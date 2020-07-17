args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  main_folder <- "."
} else {
  main_folder <- args[1]  
}

# =============================================================================
# Libraries

library(pomp)
stopifnot(packageVersion("pomp")>="2")

library(doRNG)
library(foreach)
library(doParallel)

library(plyr)
library(reshape2)
library(magrittr)

library(ggplot2)
theme_set(theme_bw())

# =============================================================================
# Parameters

output_folder <- file.path(main_folder, "output") 
if (!dir.exists(output_folder)) dir.create(output_folder)
cooking_folder <- file.path(main_folder, "output", "bake")
if (!dir.exists(cooking_folder)) dir.create(cooking_folder)
plotting_folder <- file.path(main_folder, "output", "plots")
if (!dir.exists(plotting_folder)) dir.create(plotting_folder)

input_params <- read.csv(file.path(main_folder, "input_params.csv"))

# Duration of the outbreak in units of delta_min
total_time <- input_params$T

# Total number of participants, and number of initial cases
pop_size <- input_params$N
inf0 <- input_params$I0

time_step <- 1/input_params$dT

num_sims <- 100

free_params <- c("Beta", "Gamma", "Rho")
free_param_box <- rbind(
  Beta = c(0.1, 5),
  Gamma = c(0, 3),
  Rho = c(0.5, 1)
)

fixed_params <- c("N", "S0", "I0", "R0")
fixed_param_values <- c(N = pop_size, S0 = pop_size - inf0, I0 = inf0, R0 = 0)

all_params <- c(free_params, fixed_params)

# Random seeds, keep unchanged to ensure reproducibilty of results
test_mle_seed <- 0998468235L
full_mle_seed <- 0998468235L
global_search_seed <- 290860873
test_sim_seed <- 157999
full_sim_seed <- 157999
mcap_plik_seed <- 290860873

# =============================================================================
# Observed case count

csv_table <- read.csv(file.path(main_folder, "case_counts.csv"))
case_data <- data.frame(time = csv_table$Time, cases = csv_table$Count)
max_time <- max(case_data$time)

# =============================================================================
# Csnippets defining the SIR model

# Note the use of the accumulation variable H in the measurement model. 
# C counts the number of removed new cases, not the total represented by R.
# This is the observed variable that we fit the model to.

sir_step <- Csnippet("
  double dSI = rbinom(S, 1 - exp(-Beta*I/N*dt));
  double dIR = rbinom(I, 1 - exp(-Gamma*dt));

  S -= dSI;
  I += dSI - dIR;
  R += dIR;
  C += dIR;
")

sir_init <- Csnippet("
  S = nearbyint(S0);
  I = nearbyint(I0);
  R = nearbyint(R0);
  C = 0;
")

rmeas <- Csnippet("
  // Poisson measurement model
  //double tol = 1e-6;
  //cases = rpois(Rho * C + tol);

  // Binomial measurement model
  cases = rbinom(C, Rho);

  // Overdispersed binomial model
  //double m = Rho * C;
  //double v = m * (1.0 - Rho + Psi * Psi * m);
  //double tol = 1.0e-18;
  //cases = rnorm(m, sqrt(v) + tol);
  //if (cases > 0.0) {
  //  cases = nearbyint(cases);
  //} else {
  //  cases = 0.0;
  //}
")

dmeas <- Csnippet("
  // Poisson measurement model
  //double tol = 1e-6;
  //lik = dpois(cases, Rho * C + tol, give_log);

  // Binomial measurement model
  lik = dbinom(cases, C, Rho, give_log);

  // Overdispersed binomial model
  //double m = Rho * C;
  //double v = m * (1.0 - Rho + Psi * Psi * m);
  //double tol = 1.0e-18;
  //if (cases > 0.0) {
  //  lik = pnorm(cases + 0.5, m, sqrt(v) + tol, 1, 0) - pnorm(cases - 0.5, m, sqrt(v)  + tol, 1, 0) + tol;
  //} else {
  //  lik = pnorm(cases + 0.5, m, sqrt(v) + tol, 1, 0) + tol;
  //}
  //if (give_log) lik = log(lik);
")

# =============================================================================
# POMP model

case_data %>% 
  pomp(t0 = case_data$time[1],
       time = "time",
       rprocess = euler(sir_step, delta.t=time_step),
       rinit = sir_init,
       rmeasure = rmeas,
       dmeasure = dmeas,
       accumvars = c("C"),
       statenames = c("S", "I", "R", "C"),
       partrans = parameter_trans(
         log = c("Beta", "Gamma"),
         logit = c("Rho")),
       paramnames = c(free_params, fixed_params)
  ) -> model

# =============================================================================
# IF parameters

num_test_runs <- 10

num_guesses <- 100       # Number of starting points for the parameter guesses
num_filter_iter <- 50    # Number of filtering iterations to perform
num_particles <- 2000    # Number of particles to use in filtering.
num_replicates <- 10     # Number of replicated particle filters at each point estimate

perturb_sizes <- list(Beta=0.02, Gamma=0.02, Rho=0.02)
cool_frac <- 0.5
cool_type <- "geometric"

# Variables to use in the scatterplot matrix showing the result of the IF search
pair_vars <- ~loglik+Beta+Gamma+Rho

# =============================================================================
# Test run from single starting point in parameter space and no replicates

registerDoParallel()
#set.seed(test_mle_seed, kind="L'Ecuyer")

guess <- apply(free_param_box, 1, function(x) runif(1, x[1], x[2]))

bake(file=file.path(cooking_folder, "box_search_local.rds"), {
  foreach(i=1:num_test_runs,
          .packages='pomp', .combine=c, .options.multicore=list(set.seed=TRUE)
  ) %dopar%  
  {
    mif2(
      model,
      params=c(guess, fixed_param_values),
      Np=num_particles,
      Nmif=num_filter_iter,
      cooling.type=cool_type,
      cooling.fraction.50=cool_frac,
      rw.sd=do.call(rw.sd, perturb_sizes)
    )
  }
}) -> mifs_test

ggplot(data=melt(traces(mifs_test)),
       aes(x=iteration, y=value, group=L1, color=factor(L1))) +
  geom_line() +
  guides(color=FALSE) +
  facet_wrap(~variable, scales="free_y") +
  theme_bw()

ggsave(file.path(plotting_folder, "1-mle_local_search.pdf"))

# =============================================================================
# Full MLE with multiple starting points for the free parameters

registerDoParallel()
set.seed(full_mle_seed, kind="L'Ecuyer")

stew(file=file.path(cooking_folder, "box_search_global.rda"), {
  param_guesses <- as.data.frame(apply(free_param_box, 1, function(x) runif(num_guesses, x[1], x[2])))
  
  workers <- getDoParWorkers()
  systime <- system.time({
    res_global <- foreach(guess=iter(param_guesses,"row"), 
                         .packages='pomp', .combine=rbind, .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      mifs_local <- mif2(
        model,
        params=c(unlist(guess), fixed_param_values),
        Np=num_particles,
        Nmif=num_filter_iter,
        cooling.type=cool_type,
        cooling.fraction.50=cool_frac,
        rw.sd=do.call(rw.sd, perturb_sizes)
        )
      ll <- logmeanexp(replicate(num_replicates, logLik(pfilter(mifs_local, Np=num_particles))), se = TRUE)
      data.frame(as.list(coef(mifs_local)), loglik=ll[1], loglik.se=ll[2])
    }
  })
}, seed=global_search_seed, kind="L'Ecuyer")

res_global <- as.data.frame(res_global)

all <- ldply(list(guess=param_guesses, result=subset(res_global, loglik > max(loglik)-50)), .id="type")

pairs(pair_vars, data=all, col=ifelse(all$type=="guess", grey(0.5), "red"), pch=16)
dev.copy(pdf, file.path(plotting_folder, "2-mle_global_search.pdf"))
dev.off()

# =============================================================================
# Getting best parameters

mle_params <- arrange(rbind(res_global), -loglik)
write.csv(mle_params, file=file.path(output_folder, "param_mle_global_search.csv"), row.names=FALSE, na="")

log_idx <- length(mle_params) - 1
mle_global <- mle_params[which.max( mle_params[,log_idx] ), ] 
mle_global %>% extract(all_params) %>% unlist() -> best_param_values
write.csv(best_param_values, file=file.path(output_folder, "param_point_estimates.csv"), row.names=TRUE, na="")

theta <- c(best_param_values, fixed_param_values)

# =============================================================================
# Running simulations using the MLE parameters

#set.seed(test_sim_seed)

model  %>%
  simulate(params=theta, nsim=9, format="data.frame", include.data=TRUE) -> sim_data

sim_data %>%
  ggplot(aes(x=time, y=cases, group=.id, color=(.id=="data"))) +
  guides(color=FALSE) +
  geom_line() + facet_wrap(~.id, ncol=2)

ggsave(file.path(plotting_folder, "3-simulations.pdf"))

# =============================================================================
# Comparing cumulative actual data with real data

set.seed(full_sim_seed)

model %>% 
  simulate(params=theta, nsim=num_sims, format = "data.frame", include.data=TRUE) -> sim_data

# Getting the cumulative curve for observed data
sum_data = 0
cum_data = c()
for (i in 1:max_time) {
  sum_data = sum_data + case_data$cases[i]
  cum_data = c(cum_data, sum_data)
}

# Getting the average cumulative curve for the simulations
total_size = c()
for (i in 1:num_sims) {
  temp = subset(sim_data, .id == i)
  s = sum(temp$cases)
  total_size = c(total_size, s)
}

# Taking the median
sel_sim = 0.5 * num_sims
median_index = order(total_size)[sel_sim]

sum_int = 0
cum_sim = c()
for (i in 1:max_time) {
  sim_int <- subset(sim_data, .id == median_index)
  sum_int = sum_int + sim_int$cases[i]
  cum_sim = c(cum_sim, sum_int)
}

cum_numbers = data.frame('time' = seq(1, max_time), 
                         'actual_data' = cum_data,
                         'sim_data' = cum_sim)

ggplot(cum_numbers, aes(time)) + 
  geom_line(aes(y = actual_data, colour = "Real Data")) + 
  geom_line(aes(y = sim_data, colour = "Simulated")) +
  ylab('Cumulative Number of Cases')

ggsave(file.path(plotting_folder, "4-cumulative_cases.pdf"))

# =============================================================================

# CALCULATION OF THE CONFIDENCE INTERVAL FOR THE PARAMETERS

# =============================================================================
# MCAP settings

mcap_confidence <- 0.95  # The desired confidence

# Profile design settings
mcap_profdes_len <- 15   # Number of subdivisions of the parameter range in the profile
mcap_nprof <- 10         # Number of starts per profile point
mcap_replicates <- 5
mcap_num_particles_pfilt <- 2000

# IF settings to generate the likelihood surface for each parameter
mcap_num_particles <- 1000
mcap_num_filter_iter <- 50

mcap_cool_type <- "geometric"
mcap_cool_frac <- 0.5
mcap_cool_frac_lastif <- 0.1

mcap_lambda <- c(Beta=0.75, Gamma=0.75, Rho=0.75) # Lambda closer to 1 increases the curvature of the likelihood approximation
mcap_ngrid <- c(Beta=1000, Gamma=1000, Rho=1000)

# =============================================================================
# Function to calculate the MCAP CIs 

mcap <- function(lp, parameter, confidence, lambda, Ngrid) {
  smooth_fit <- loess(lp ~ parameter, span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit, newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter - smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weight <- rep(0, length(parameter))
  weight[included] <- (1 - (dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weight=weight,
                      data = data.frame(lp=lp, b=parameter, a=-parameter^2))
  b <- unname(coef(quadratic_fit)["b"] )
  a <- unname(coef(quadratic_fit)["a"] )
  m <- vcov(quadratic_fit)
  
  var_b <- m["b", "b"]
  var_a <- m["a", "a"]
  cov_ab <- m["a", "b"]
  
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1/(2*a)
  se_total_squared <- se_mc_squared + se_stat_squared
  
  delta <- qchisq(confidence, df=1) * ( a * se_mc_squared + 0.5)
  loglik_diff <- max(smoothed_loglik) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(lp=lp, parameter=parameter, confidence=confidence,
       quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
       smooth_fit=smooth_fit,
       fit=data.frame(parameter=parameter_grid,
                      smoothed=smoothed_loglik,
                      quadratic=predict(quadratic_fit, list(b = parameter_grid, 
                                                            a = -parameter_grid^2))),
       mle=smooth_arg_max, ci=ci, delta=delta,
       se_stat=sqrt(se_stat_squared), se_mc=sqrt(se_mc_squared), 
       se=sqrt(se_total_squared)
  )
}

# =============================================================================
# Generates the model for the profile likelihood calculation

plik <- function(pname, pval0, pval1) {
  registerDoParallel()
  
  bake(file=file.path(cooking_folder, paste(pname, "_mcap.rds", sep="")), {  
    
    desel <- which(names(mle_params) %in% c("loglik", "loglik.se", pname))
    mle_params %>% 
      subset(
        loglik > max(loglik,na.rm=TRUE) - 20,
        select=-desel
      ) %>% 
      melt(id=NULL) %>% 
      daply(~variable, function(x)range(x$value)) -> box
    
    starts <- profileDesign(pname=seq(pval0, pval1, length=mcap_profdes_len),
                            lower=box[,1], upper=box[,2],
                            nprof=mcap_nprof)
    names(starts) <- sub("pname", pname,names(starts))
    
    psizes <- perturb_sizes[names(perturb_sizes) != pname]
    
    foreach(params=iter(starts, "row"),
            .combine=rbind,
            .packages="pomp",
            .options.multicore=list(set.seed=TRUE),
            .options.mpi=list(seed=mcap_plik_seed, chunkSize=1)
    ) %dopar% {
      mf <- mif2(model,
                 params=unlist(params),
                 Np=mcap_num_particles,
                 Nmif=mcap_num_filter_iter,
                 cooling.type=mcap_cool_type,
                 cooling.fraction.50=mcap_cool_frac,
                 rw.sd=do.call(rw.sd, psizes)
      )
      mf <- mif2(mf, 
                 Np=mcap_num_particles,
                 Nmif=mcap_num_filter_iter,
                 cooling.fraction.50=mcap_cool_frac_lastif)
      ll <- logmeanexp(replicate(mcap_replicates, 
                                 logLik(pfilter(mf, Np=mcap_num_particles_pfilt))), se=TRUE)
      data.frame(as.list(coef(mf)), loglik=ll[1], loglik.se=ll[2])
    }
  }) -> pmodel
  
  return(pmodel)
}

# =============================================================================
# Iterate over all the free parameters in the model to calculate their CIs... may take a while!

par_names <- row.names(free_param_box)
for (i in 1:nrow(free_param_box)) {
  name <- par_names[i]  
  
  row <- free_param_box[i,]
  mdl <- plik(pname=name, pval0=row[1], pval1=row[2])
  
  par_range <- seq(row[1], row[2], length=mcap_profdes_len)
  log_likelihoods <- c()
  for (val in par_range) {
    likelihoods = subset(mdl, abs(mdl[[name]]-val)<1)$loglik
    if (length(likelihoods)==0) next
    log_likelihoods = c(log_likelihoods, max(likelihoods))
  }
  
  x <- mcap(log_likelihoods, par_range, mcap_confidence, mcap_lambda[[i]], mcap_ngrid[[i]])
  if (i == 1) {
    cis <- data.frame("name" = c(name), "x0" = c(x$ci[1]), "x1" = c(x$ci[2]), stringsAsFactors = FALSE)  
  } else {
    cis <- rbind(cis, c(name, x$ci[1], x$ci[2]))  
  }

  ggplot(x$fit, aes(parameter, quadratic)) + geom_line() + 
    geom_vline(xintercept=c(x$ci[1], x$ci[2]), linetype=4, colour='red') +
    geom_point(data = data.frame('parameters'=par_range, 'loglik'=log_likelihoods), 
               aes(parameters, log_likelihoods)) 
  ggsave(file.path(plotting_folder, paste("5-", name, "_ci.pdf", sep="")))
}

write.csv(cis, file=file.path(output_folder, "param_confidence_intervals.csv"), row.names=FALSE, na="")
