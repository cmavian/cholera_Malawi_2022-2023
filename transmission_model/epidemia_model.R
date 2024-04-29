library(epidemia)
library(EpiEstim)
library(lubridate)
library(rstanarm)
library(tidyverse)
library(patchwork)
library(harrypotter)
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
rt_plot  <- spaghetti_rt(fm1) + ylim(c(0.2,1.6)) +
  ylab(TeX("$\\R_t$ estimate")) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  theme(axis.text.x=element_text(angle=0, hjust=1))

# Plot Observed vs Estimates cases
obs_plot <- spaghetti_obs(fm1,type="Cases") + ylab("Observed and Fitted Cases") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  theme(axis.text.x=element_text(angle=0, hjust=1)) + theme(legend.position = "none")

# Plot effect sizes for covariates
pars_plot <- plot(fm1, par_models = "R", par_types = "fixed") # ,prob_outer = 0.8


(rt_plot + obs_plot)
pars_plot

pp_40 <- pars_plot

if(0)
{
  rt_plot_fv <- rt_plot
  obs_plot_fv <- obs_plot
  pars_plot_fv <- pars_plot
}

if(SAVE_PLOTS)
{
  ggsave("rt_estimate.pdf",rt_plot,device = cairo_pdf)
  ggsave("cases_plot.pdf",obs_plot,device = cairo_pdf)
  ggsave("effect_sizes.pdf",pars_plot,device = cairo_pdf)
}


#overview plot
#overview_plot <- (rt_plot_fv + obs_plot_fv + pars_plot_fv) / (rt_plot_f + obs_plot_f + pars_plot_f) / (rt_plot_no_covariates + obs_plot_no_covariates + pars_plot_no_covariates)
#ggsave("overview_plot.pdf",overview_plot,device = cairo_pdf, width=20, height=16)

# extract Rt estimates from epidemia
Rt_adj <- colMeans(posterior_rt(fm1)$draw)

data$Rt_adj <- Rt_adj
data %>% dplyr::select(date,Rt_adj,Cases,ma_cases,rain_proxy,vaccination) %>% mutate('Cases/1000'=Cases/1000,ma_cases=ma_cases/1000) %>%
  dplyr::select(-Cases) %>%
  pivot_longer(-date,names_to = "time_series",values_to = 'value') %>%
  group_by(time_series) %>% ggplot(aes(x=date,y=value,col=time_series)) + geom_line() +
  theme_light() + scale_color_hp_d(option="NewtScamander") + ylab("")

#new_data <- get_cholera_data() %>% filter(date > date[which(cumsum(Deaths) > 10)[1] - 30])
#plot_obs(fm1, type = "Cases", newdata = newdata, groups = "Malawi")

# We can validate our estimates using EpiEstim
res_parametric_si <- estimate_R(data.frame(dates=data$date,I=data$Cases),
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 5,
                                  std_si = 8))
)

epiestim_compare <- plot(res_parametric_si, 'R') + geom_line(y=Rt_adj[8:length(Rt_adj)],col='blue')

if(SAVE_PLOTS)
{
  ggsave("EpiEstim_epidemia_compare.pdf",epiestim_compare,device = cairo_pdf)
}

