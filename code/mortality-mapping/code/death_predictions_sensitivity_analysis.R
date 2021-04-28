library(data.table)   
library(ggplot2)
library(dplyr)
library(cowplot)

# set the directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# pull in functions for death predictions
source('code/covid-19-data-functions.R')

make_plot <- function(res, title, china_br = F, china_mr = F, scaled_true = T, legend = F, 
                      min_date = as.Date('2020-02-01')){
  # only use past the min date
  res = res[res$date >= min_date,]
  
  a = data.frame(date = res$date, cases = res$cum_cases_estimate, type = 'mMAP',
                 stringsAsFactors = F)
  b = data.frame(date = res$date, cases = cumsum(res$daily_cases_true), type = 'reported',
                 stringsAsFactors = F)
  
  # if adding in the true scaled plots
  if(scaled_true){
    c = b
    c$cases = b$cases*mean(a$cases)/mean(b$cases)
    c$type = 'reported-scaled'
    
    tt = rbind(a, b, c)
  }else{
    tt = rbind(a, b)
  }
  
  
  if(legend){
    p1 <- ggplot(data = tt, mapping = aes(x = date, y = cases, col = type)) +
      geom_line(size = 1) +
      scale_color_manual(values = c('#a167c9','black','grey60')) +
      ggtitle(title) + 
      theme(legend.box.margin = margin(0, 0, 0, 12),
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.position = 'bottom',
            legend.text = element_text(size = 12)) +
      guides(color = guide_legend(override.aes = list(size=6)))
  }else{
    p1 <- ggplot(data = tt, mapping = aes(x = date, y = cases, col = type)) +
      geom_line(size = 1) +
      scale_color_manual(values = c('#a167c9','black','grey60')) + 
      ggtitle(title) +
      theme_minimal() +
      theme(legend.position = 'none',
            axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + 
      aes(group=rev(type))
    if(title == 'China'){
      if(china_br){
        p1 = p1 + theme(plot.title = element_text(vjust = -20, hjust = 0.8))
      }else if(china_mr){
        p1 = p1 + theme(plot.title = element_text(vjust = -5, hjust = 0.05))
      }
    }else{
      p1 = p1 + theme(plot.title = element_text(vjust = -10, hjust = 0.1))
    }
  }
  
  return(p1)
}

countries = c('US','Italy','Spain','Germany', 'Japan', 'Korea, South')

#
##### (1) Simulate death counts and then back-predict true case counts #####
#cases = get_death_predicted_case_data(level = 'US', min_deaths = 20)
# get deaths by country
deaths_reported = get_time_series_data_JHU(type = 'deaths_global', US_only = F)

# get cases
cases_list = get_time_series_data_JHU(type = 'confirmed_global', US_only = F)
cases = cases_list[[1]]

# date to right censor at
# right_censor_date = as.Date(cases_list[[2]], format = '%m-%d-%Y') + 50
# right_censor_date = as.Date(cases_list[[2]], format = '%m-%d-%Y')
right_censor_date = as.Date('2020-06-07')

# only keep countries with > 1000 cases

# get probabilities of deaths
death_df = return_onset_death_df()
p_death = death_df$p

# infected fatality rate
IFR = 0.013

# initialize plot list
plot_list = list()
iter = 1

set.seed(1)
warning('simulation is random - I could introduce an average over runs')
# cycle through set of countries
for(country in countries){
  print(country)
  # get the cases for that country
  C = get_case_data(cases, state = NULL, country = country)
  
  # simulate deaths going forward
  deaths = simulate_deaths(cases = C, p_death = p_death, IFR = IFR)
  
  # right censor
  deaths = deaths[as.Date(names(deaths)) <= right_censor_date]
  deaths_2 = data.frame(date = as.Date(names(deaths)), death = deaths)
  
  # back-track deaths to cases
  res = deaths_to_cases(deaths_2, p_death, IFR = IFR, stopping_criteria = 'chi-2')$cases
  
  # get cumulative cases
  res$cum_cases_estimate = cumsum(res$daily_cases_estimate)
  
  # cut off the last week
  res = res[res$date <= right_censor_date - 7,]
  
  # add in original cases to data frame:
  res$daily_cases_true = 0
  C = C[as.Date(names(C)) %in% res$date]
  res$daily_cases_true[match(as.Date(names(C)), res$date)] = C
  
  if(country == 'Korea, South'){country = 'South Korea'}
  if(country == 'US'){country = 'United States'}
  
  plot_list[[iter]] = make_plot(res, title = country, china_br = T, scaled_true = F)
  iter = iter + 1
}

legend = get_legend(make_plot(res, country, china_br = T, scaled_true = F, legend = T))
plot_grid(plot_grid(plotlist = plot_list),  legend, nrow = 2, ncol = 1, rel_heights = c(3,.4))
# ggsave('results/final_sim_death_plots_11132020.png', width = 10, height = 5)

# read in cases from random countries
# repeatedly sample potential death distributions given IFR and death distr
# using death distribution, back-predict cases
  # and the current date cutoff for right censoring!
# compare cases vs. uncovered cases - MSE, correlation, mean overall


##### (2) Compare predicted cases to reported cases #####
#cases = get_death_predicted_case_data(level = 'US', min_deaths = 20)
# get deaths by country
deaths_list = get_time_series_data_JHU(type = 'deaths_global', US_only = F)
death_data = deaths_list[[1]]
cases_list = get_time_series_data_JHU(type = 'confirmed_global', US_only = F)
cases = cases_list[[1]]

# # date to right censor at
# right_censor_date = as.Date(cases_list[[2]], format = '%m-%d-%Y')
# right_censor_date = as.Date('2020-06-07')

# countries with good testing
# countries = c('Germany', 'Korea, South') #, 'Czechia', 'Norway', 'Japan') # , 'Taiwan*', 'Singapore') <- too few deaths

warning('removing China')
countries = c('US','Italy','Spain','Germany', 'Japan', 'Korea, South')
#countries = c('United Kingdom', 'France', 'Netherlands', 'Japan', 'Taiwan')


# get probabilities of deaths
death_df = return_onset_death_df()
p_death = death_df$p

# infected fatality rate
IFR = 0.013

# smoothing distance for deaths
smoothing_dist = 3

# initialize
plot_list = list()
iter = 1

# cycle through set of countries
for(country in countries){
  print(country)
  # get the cases for that country
  C = get_case_data(cases, state = NULL, country = country)
  
  # get the true deaths reported
  deaths = get_case_data(death_data, state = NULL, country = country, smoothing_dist = smoothing_dist)
  print(sum(deaths))
  
  # right censor
  deaths = deaths[as.Date(names(deaths)) <= right_censor_date]
  
  deaths_2 = data.frame(date = as.Date(names(deaths)), death = deaths)
  
  # back-track deaths to cases
  res = deaths_to_cases(deaths_2, p_death, IFR = IFR, stopping_criteria = 'chi-2')$cases
  
  # get cumulative cases
  res$cum_cases_estimate = cumsum(res$daily_cases_estimate)
  
  # cut off the last week
  res = res[res$date <= right_censor_date - 7,]
  
  # add in original cases to data frame:
  res$daily_cases_true = 0
  C = C[as.Date(names(C)) %in% res$date]
  res$daily_cases_true[match(as.Date(names(C)), res$date)] = C
  
  if(country == 'Korea, South'){country = 'South Korea'}
  if(country == 'US'){country = 'United States'}
  
  plot_list[[iter]] = make_plot(res, country, china_mr = T, legend = F)#(iter %% 3 == 0))
  iter = iter + 1
}

legend = get_legend(make_plot(res, country, legend = T))
plot_grid(plot_grid(plotlist = plot_list),  legend, nrow = 2, ncol = 1, rel_heights = c(3,.4))

# ggsave('results/final_true_death_plots_11132020.png', width = 10, height = 5)
