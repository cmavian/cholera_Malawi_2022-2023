### Run Transmission model ###

if(!require(epidemia))         # Check that epidemia is installed, if not follow below link
{
  errorCondition('epidemia installation is required\n 
                 #install.packages("devtools")
                 devtools::install_github("ImperialCollegeLondon/epidemia")')
}

if(!require(rstan)|!require(rstanarm))    # Check that rstan and rstanarm are installed, if not follow below link
{
  errorCondition('rstan installation is required, please these instructions:\n
                 https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started')
}

library(epidemia)
library(EpiEstim)
library(lubridate)
library(rstanarm)
library(tidyverse)
library(patchwork)
library(zoo)
library(latex2exp)

source('transmission_model/helper_functions.R')

options(mc.cores = parallel::detectCores())
SAVE_PLOTS <- FALSE

# Get Data
data <- get_cholera_data()
# We determine the start of the epidemic (for estimation of the transmission model) as 30 days
# before the the 5th cumulative death
data <- data %>%
  filter(date > date[which(cumsum(Deaths) > 5)[1] - 30])

# Generate the Serial Interval and Infection to Onset distributions
cholera_SI  <- get_cholera_SI()
cholera_i2o <- get_cholera_i2o()

# Transmission: define linear predictors and random walk
# weekly random walk to obtain more stable estimates
rt <- epirt(formula = R(Country, date) ~ 1 + rain_proxy + vaccination + rw(prior_scale = 0.01,time=week),
            prior = rstanarm::normal(scale = 1),
            prior_intercept = normal(log(1.56), 0.2), link = 'log')

# Observations: define observations (in our case we only use cases) and the infection to onset distribution
obs <-  epiobs(formula = Cases ~ 1, link = "identity",
               i2o = cholera_i2o/sum(cholera_i2o))

# Infections: define renewal equation to propagate infections
# define initial population size, adjust population for cases and define prior for susceptible population
inf     <- epiinf(gen = cholera_SI/sum(get_cholera_SI()), pop_adjust = TRUE, pops = pop,
                  prior_susc = normal(.95, 0.1))

args <- list(rt = rt, obs = obs, inf = inf, data = data, iter = 2e3, seed = 12345)
args_ext <- args

# Estimate epidemia model
system.time(fm1 <- do.call(epim, args))

# Plot Rt estimates 
## Figure 3(B)    inset above case curve
## Figure S3(A)   [if model is run without covariates this would be S3(C)]
rt_plot  <- spaghetti_rt(fm1) + ylim(c(0.2,1.6)) +
  ylab(TeX("$\\R_t$ estimate")) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  theme(axis.text.x=element_text(angle=0, hjust=1))

# Plot Observed vs Estimates cases
## Figure S3(B)
obs_plot <- spaghetti_obs(fm1,type="Cases") + ylab("Observed and Fitted Cases") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  theme(axis.text.x=element_text(angle=0, hjust=1)) + theme(legend.position = "none")

# Plot effect sizes for covariates
## Figure S3(D) [generate for different priors of S0 to reproduce full subplot D]
pars_plot <- plot(fm1, par_models = "R", par_types = "fixed") 


if(SAVE_PLOTS)
{
  ggsave("rt_estimate.pdf",rt_plot,device = cairo_pdf)
  ggsave("cases_plot.pdf",obs_plot,device = cairo_pdf)
  ggsave("effect_sizes.pdf",pars_plot,device = cairo_pdf)
}

CHECK_MODEL_WITH_EPIDEMIA = FALSE
if(CHECK_MODEL_WITH_EPIDEMIA)
{
  # extract Rt estimates from epidemia
  Rt_adj <- colMeans(posterior_rt(fm1)$draw)
  
  # We can validate our estimates using EpiEstim Rt estimation
  res_parametric_si <- estimate_R(data.frame(dates=data$date,I=data$Cases),
                                  method="parametric_si",
                                  config = make_config(list(mean_si = 5, std_si = 8))
  )
  
  # plot comparison of EpiEstim and epidemia Rt estimation
  epiestim_compare <- plot(res_parametric_si, 'R') + geom_line(y=Rt_adj[8:length(Rt_adj)],col='blue')
}
