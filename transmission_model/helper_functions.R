### Helper functions to estimate transmission model

#' Infection to onset distribution
#'
#' @return i2o distribution
get_cholera_i2o <- function()
{
  set.seed(1234) # set for reproducibility

  x <- seq(0.1, 10 , by=0.1)
  f <- as_tibble(list(x=x-0.1,p=dlnorm(x,log(1.4),log(1.98))))
  f$floor <- floor(f$x)
  out <- f %>% group_by(floor) %>% summarize(prob=sum(p*0.1))

  return(out$prob)
}

#' Generate Serial Interval
#'
#' @param plot_SI    Plot pdf / cdf of serial interval for validation, default to FALSE
#' @param alpha      alpha parameter of Gamma distribution
#' @param theta      theta parameter of Gamma distribution
#' @param n          number of points to generate empirical distribution
#' @param shift      how many days to shift
#'
#' @return Serial Interval
get_cholera_SI <- function(plot_SI=FALSE,alpha=3.125,theta = c(0.625),n=100,shift=1)
{
  set.seed(1234) # set for reproducibility

  n = n - shift

  x <- seq(0, n-1 , len=n)

  f <- sapply(theta, function(mytheta) dgamma(x, shape = alpha, rate = mytheta))
  F <- sapply(theta, function(mytheta) pgamma(x, shape = alpha, rate = mytheta))
  leg <- sapply(theta, function(mytheta) sprintf("G(%.3g, %.3g)", alpha, mytheta))

  if(plot_SI)
  {
    par(mfrow=c(1,2))

    kk <- seq(dim(f)[2]); cpal <- rainbow(max(kk))

    plot(x, f[,1], col=cpal[1],
         type="l", lwd=1,
         xlim=c(0.9*min(x), 1.1*max(x)), ylim=c(0, 1.1*max(f)),
         xlab="x", ylab="f(x)", main="PDF of Gamma distributions")

    for(k in kk[-1]) lines(x, f[,k], col=cpal[k], lwd=1)
    legend("topright", col=cpal, lwd=1, legend=leg)


    plot(x, F[,1], col=cpal[1], type="l", lwd=1,
         xlim=c(0.9*min(x), 1.1*max(x)), ylim=c(0, 1.1*max(F)),
         xlab="x", ylab="F(x)", main="CDF of continuous random variables")

    for(k in kk[-1]) lines(x, F[,k], col=cpal[k], lwd=1)
    legend("bottomright", col=cpal, lwd=1, legend=leg)
  }

  return(c(rep(0,shift),f[,1]))
}

#' Generate tibble with all required data for transmission model
#'
#' @return data tibble with cases, deaths, flooding, vaccination and population data
get_cholera_data <- function()
{
  malawi_cases       <- read_csv("transmission_model/cases.csv")
  malawi_deaths      <- read_csv("transmission_model/deaths.csv") %>% rename('deaths'='Cases')

  flooding      <- read_csv("transmission_model/flooding_malawi.csv") %>% mutate(day_end=as.Date(as.character(day_end),"%Y%m%d"))

  cholera_data  <- malawi_cases %>% left_join(malawi_deaths,by=c("Date"))

  date <- seq(as.Date(min(cholera_data$Date)), as.Date(max(cholera_data$Date)), by="days")
  data <- data.frame(Country = "Malawi", fill_cases = rep(NA,length(date)),
                     Date = date)

  cholera_data <- data %>% left_join(cholera_data,by=c('Date')) %>%
    left_join(flooding %>% dplyr::select("day_end","both_layers_modified"),by=c("Date"="day_end")) %>%
    rename(flooding=both_layers_modified) %>% fill(flooding) %>% dplyr::select(-fill_cases) %>%
    rename(Deaths=deaths,date=Date) %>%
    mutate(Cases            = replace_na(Cases,0),
           Deaths           = replace_na(Deaths,0),
           flooding         = flooding/10000,
           ma_cases         = rollapply(Cases,7,mean,align='right',fill=NA),
           sign_diff_cases  = sign(ma_cases-dplyr::lag(ma_cases,7)),
           diff_flooding    = replace_na(flooding-dplyr::lag(flooding,1),0),
           rain_proxy       = (dplyr::lag(diff_flooding * (diff_flooding>0),7)>0.1),  # 5 + 1.4 = mean SI + incubation period
           diff_flooding_l2 = replace_na(flooding-dplyr::lag(flooding,2),0),
           diff_flooding_l3 = replace_na(flooding-dplyr::lag(flooding,3),0),
           diff_flooding_l4 = replace_na(flooding-dplyr::lag(flooding,4),0),
           weekday          = lubridate::wday(date,label=TRUE))

  # need to add vaccination
  cholera_data$vaccination <- 0
  cholera_data[cholera_data$date>as.Date("2022-10-31") & cholera_data$date < as.Date("2023-02-01"), ]$vaccination <- seq(0,1,by=1/91) #60
  cholera_data[cholera_data$date >= as.Date("2023-02-01"), ]$vaccination <- 1
  cholera_data <- cholera_data %>% mutate(vaccination=dplyr::lag(vaccination,10)) # 10 days for vaccination to be effective

  cholera_data$week    <- floor(seq(1:dim(cholera_data)[1])/7)
  cholera_data$pop     <- 20091635

  return(as_tibble(cholera_data))
}
