rm(list=ls())
library(data.table)   
library(ggplot2)
library(dplyr)
library(cowplot)

# set the directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# pull in functions for death predictions
source('code/covid-19-data-functions.R')

##### Getting 4 different estimates for mMAP - 11/10/2020 #####

# OG unadjusted (not used in paper)
D_unadj_old = get_death_predicted_case_data(level = 'US', min_deaths = 0, asymptomatic_scale_down = 0.4, IFR_scale = 1,  add_pneumonia = F, remove_last = 5, stopping_criteria = 'chi-2')

# new unadjusted
D_unadj = get_death_predicted_case_data(level = 'US', min_deaths = 0, asymptomatic_scale_down = 0.5, IFR_scale = 1,  add_pneumonia = F, remove_last = 5, stopping_criteria = 'chi-2')

# adjusted - no pneumonia (not used in paper)
D_adj_noPneum = get_death_predicted_case_data(level = 'US', min_deaths = 0, asymptomatic_scale_down = 0.4, IFR_scale = 0.65/1.1,  add_pneumonia = F, remove_last = 5, stopping_criteria = 'chi-2')

# adjusted - pneumonia
D_adj = get_death_predicted_case_data(level = 'US', min_deaths = 0, asymptomatic_scale_down = 0.4, IFR_scale = 0.65/1.1,  add_pneumonia = T, remove_last = 5, stopping_criteria = 'chi-2')

# create weekly cumulative estimates

setdiff(D_unadj_old$df$location, D_unadj$df$location)
setdiff(D_unadj_old$df$location, D_adj$df$location)
length(unique(D_unadj_old$df$location))
length(unique(D_unadj$df$location))
length(unique(D_adj$df$location))
length(unique(D_adj_noPneum$df$location))
# ok good.

# convert to weekly data for cases and prediction
# location | week-end | N | method

y = 'cum_cases_estimate'
df = NULL
weekstarts = as.Date(c('2020-03-01')) + seq(0, length.out = 11, by =7)
weekendings = as.Date(c('2020-03-01')) + seq(6, length.out = 11, by =7)

# cycle through each method
for(method in c('D_unadj_old', 'D_unadj', 'D_adj', 'D_adj_noPneum')){
  death_pred = get(method)
  
  for(loc in unique(death_pred$df$location)){
    tmp = death_pred$df %>% filter(location == loc)
    
    # tmp data frame to store values
    tmp_2 = data.frame(location = loc,
                       week_ending = weekendings,
                       week_start = weekstarts,
                       N = NA, 
                       method = method,
                       stringsAsFactors = F)
    
    # get the values for the week
    for(i in 1:nrow(tmp_2)){
      tmp_2[i,'N'] = tmp %>% filter(date == weekendings[i]) %>% pull(y) - 
        tmp %>% filter(date == weekendings[i]-7) %>% pull(y)
    }
    
    df = rbind(df, tmp_2)
  }
}

df %>% filter(location == 'United States') %>% group_by(method) %>% summarize(sum(N))

  
# Looking at the fits
par(mfrow=c(3,3))
for(loc in unique(D_unadj[[1]]$location)){
  tmp = D_unadj[[1]] %>% filter(location == loc)
  tmp2 = D_adj[[1]] %>% filter(location == loc)
  plot(tmp$date, tmp$daily_cases_estimate, type = 'l', main = loc)
  lines(tmp2$date, tmp2$daily_cases_estimate, type = 'l', col = 'red')
  print(loc)
  # print(sprintf('chi = %s and %s', death_pred[[2]][loc], death_pred_v2[[2]][loc]))
  print(tmp$cum_cases_estimate[nrow(tmp)-5])
  print(tmp2$cum_cases_estimate[nrow(tmp2)-5])
  # print('----')
}
# ok. Still US, NY, NJ, and TX a concern. I will look at using different data dates.

# creating final results
df2 = df %>% filter(method %in% c('D_unadj','D_adj'))
df2$method = gsub('D','mMAP',df2$method)
write.table(df2, file = 'results/mMAP_results_11132020.txt', sep=',', row.names = F, quote = F)

##### Comparing chi-squared fixed and original chi-squared results #####

D = read.table('results/mMAP_results_11122020.txt', sep = ',', header = T)
D2 = read.table('results/mMAP_results_11132020.txt', sep = ',', header = T)

D = D %>% filter(method == 'mMAP_adj')
D2 = D2 %>% filter(method == 'mMAP_adj')

# want to look at US, NY, NJ, TX, GA
for(loc in c('United States', 'New York', 'New Jersey', 'Texas', 'Georgia')){
  print(loc)
  tmp = D %>% filter(location == loc)
  tmp2 = D2 %>% filter(location == loc)
  # print(head(tmp))
  # print('---above: OG, below = NEW ----')
  # print(head(tmp2))
  print(sprintf('sum OG: %s, sum NEW: %s', sum(tmp$N), sum(tmp2$N)))
}

##### Testing national sum vs. state-level sums #####

names(D_adj)

deaths = unlist(D_adj$death_list)

# US
deaths[53]
sum(deaths[1:52])
# so a little higher

# made the df above for just D_adj
df = df %>% filter(week_ending <= '2020-04-04')

sum(df$N[df$location == 'United States'])
sum(df$N[df$location != 'United States'])

# getting deaths df
D_list = get_time_series_data_NYT(type = 'states')
D = D_list[[1]]; D_death = D_list[[2]]; date = D_list[[3]]

sum(D_death$`04/25/20`) # 106,654
sum(D_death$`04/25/20`[D_death$state == 'United States']) #53,315
sum(D_death$`04/25/20`[!(D_death$state %in% c('United States', 'Northern Mariana Islands', 'Virgin Islands','Guam'))]) #53,327

# barely any different.

# hm ok. We would expect the US counts and all other counts to converge right?

# cases
y = 'cum_cases'
df = NULL
weekstarts = as.Date(c('2020-03-01')) + seq(0, length.out = 11, by =7)
weekendings = as.Date(c('2020-03-01')) + seq(6, length.out = 11, by =7)

# cycle through each method
method  = 'D_unadj_old'
death_pred = get(method)

for(loc in unique(death_pred$df$location)){
  tmp = death_pred$df %>% filter(location == loc)
  
  # tmp data frame to store values
  tmp_2 = data.frame(location = loc,
                     week_ending = weekendings,
                     week_start = weekstarts,
                     N = NA, 
                     method = method,
                     stringsAsFactors = F)
  
  # get the values for the week
  for(i in 1:nrow(tmp_2)){
    tmp_2[i,'N'] = tmp %>% filter(date == weekendings[i]) %>% pull(y) - 
      tmp %>% filter(date == weekendings[i]-7) %>% pull(y)
  }
  
  df = rbind(df, tmp_2)
}


sum(df$N[df$location == 'United States'])
sum(df$N[df$location != 'United States'])

df2 = df %>% filter(week_start <= '2020-04-01')
sum(df2$N[df2$location == 'United States'])
sum(df2$N[df2$location != 'United States'])

# save it
df = df %>% select(-method)
write.csv(df, file = 'results/reported_cases_by_week_11182020.csv', row.names = F)

##### Testing training with different cutoff dates - 11/10/2020 #####
# new unadjusted
D_unadj = get_death_predicted_case_data(level = 'US', min_deaths = 0, asymptomatic_scale_down = 0.5, IFR_scale = 1,  add_pneumonia = F, remove_last = 5, cutoff_date = '2020-06-15')

for(loc in c('Texas', 'New Jersey', 'New York', 'United States')){
  tmp = D_unadj[[1]] %>% filter(location == loc)

  plot(tmp$date, tmp$daily_cases_estimate, type = 'l', main = loc)

  # print('----')
}

# NJ and TX converged. NY did better. US did better

# just testing
y = 'cum_cases_estimate'
df_backup = df
df = NULL
weekstarts = as.Date(c('2020-03-01')) + seq(0, length.out = 11, by =7)
weekendings = as.Date(c('2020-03-01')) + seq(6, length.out = 11, by =7)

# cycle through each method
for(method in c('D_unadj')){
  death_pred = get(method)
  
  for(loc in unique(death_pred$df$location)){
    tmp = death_pred$df %>% filter(location == loc)
    
    # tmp data frame to store values
    tmp_2 = data.frame(location = loc,
                       week_ending = weekendings,
                       week_start = weekstarts,
                       N = NA, 
                       method = method,
                       stringsAsFactors = F)
    
    # get the values for the week
    for(i in 1:nrow(tmp_2)){
      tmp_2[i,'N'] = tmp %>% filter(date == weekendings[i]) %>% pull(y) - 
        tmp %>% filter(date == weekendings[i]-7) %>% pull(y)
    }
    
    df = rbind(df, tmp_2)
  }
}

df %>% filter(location == 'United States') %>% group_by(method) %>% summarize(sum(N))
df_backup %>% filter(location == 'United States') %>% group_by(method) %>% summarize(sum(N))


df %>% filter(location == 'Texas') %>% group_by(method) %>% summarize(sum(N))
df_backup %>% filter(location == 'Texas') %>% group_by(method) %>% summarize(sum(N))

# so not that different anyway. Like 1% different. Not a big deal.


##### Testing chi-squared convergence #####
system.time({death_pred = get_death_predicted_case_data(level = 'US', min_deaths = 0, asymptomatic_scale_down = 0.4, stopping_criteria = 'chi', gonzalo_mMAP = T, remove_last = 5)})

system.time({death_pred_v2 = get_death_predicted_case_data(level = 'US', min_deaths = 0, asymptomatic_scale_down = 0.4, stopping_criteria = 'chi', gonzalo_mMAP = T, smoothing_dist = 3, remove_last = 5)})
# So this ^ almost always converges (weirdly not for the national level - I don't know what's up with that.)
# Well the first one also doesn't work for vermont currently


# Thoughts - the chi squared is measuring how close the predicted deaths from cases would be compared to reported. The reason (or one of them) why chi squared isnt reaching one is probably that the deaths are just not recoverable. They are not reasonable if deaths were actually reported consistently. It probably has a lot to do with the lack of reporting on weekends. I will try this with weekly smoothing as well and see what happens.

par(mfrow=c(3,3))
for(loc in unique(death_pred[[1]]$location)){
  tmp = death_pred[[1]] %>% filter(location == loc)
  tmp2 = death_pred_v2[[1]] %>% filter(location == loc)
  plot(tmp$date, tmp$daily_cases_estimate, type = 'l', main = loc)
  lines(tmp2$date, tmp2$daily_cases_estimate, type = 'l', col = 'red')
  print(loc)
  print(sprintf('chi = %s and %s', death_pred[[2]][loc], death_pred_v2[[2]][loc]))
  print(tmp$cum_cases_estimate[nrow(tmp)-5])
  print(tmp2$cum_cases_estimate[nrow(tmp2)-5])
  print('----')
}

# initialize death estimates
case_estimates = data.frame(location = unique(death_pred[[1]]$location), mMAP_est = 0, mMAP_smoothed = 0, simple_est = 0, mMAP_converge = NA, mMAP_smooth_converge = NA,
                            stringsAsFactors = F)

warning('no Vermont!')

# get the case estimates from each method
for(i in 1:nrow(case_estimates)){
  loc = case_estimates$location[i]
  tmp = death_pred[[1]] %>% filter(location == loc)
  tmp2 = death_pred_v2[[1]] %>% filter(location == loc)
  
  case_estimates$mMAP_est[i] = tmp$cum_cases_estimate[nrow(tmp)-5]
  case_estimates$mMAP_smoothed[i] = tmp2$cum_cases_estimate[nrow(tmp2)-5]
  # estimating cases as deaths/(IFR*asymptomatic scale down)
  case_estimates$simple_est[i] = death_pred$death_list[[loc]]/death_pred$IFR_list[[loc]]*0.6
  case_estimates$mMAP_converge[i] = (death_pred$chi_list[[loc]] <= 1)
  case_estimates$mMAP_smooth_converge[i] = (death_pred_v2$chi_list[[loc]] <= 1)
}

case_estimates$ratio = case_estimates$mMAP_smoothed/case_estimates$simple_est
case_estimates = case_estimates[order(case_estimates$ratio),]
fix(case_estimates)
# Texas, New York, New Jersey, and the United States did not converge for the smoothed method. United States is a big issue. I need to revisit that. But right now I need a probability break
save(case_estimates, death_pred, death_pred_v2, file = 'results/tmp_state_estimates_09272020.RData')

# so I'm definitely concerned about Texas and the US, but so be it.



##### Creating the final results #####
# Want Location | mMAP_Symptomatic_Cases_9_19_2020 | mMAP_Infections_9_19_2020
cases = data.frame(location = unique(death_pred_v2[[1]]$location), mMAP_Symptomatic_Cases_9_19_2020 = 0, mMAP_Infections_9_19_2020 = 0,  stringsAsFactors = F)

for(i in 1:nrow(cases)){
  loc = cases$location[i]
  tmp = death_pred_v2[[1]] %>% filter(location == loc)

  cases$mMAP_Symptomatic_Cases_9_19_2020[i] = tmp$cum_cases_estimate[nrow(tmp)-5]
  cases$mMAP_Infections_9_19_2020[i] = cases$mMAP_Symptomatic_Cases_9_19_2020[i]/0.6
}

write.csv(cases, row.names=  F, file = 'results/mMAP_case_estimates_09192020.csv')

##### Creating the boxplot #####
# read in
df = fread('results/compiled_estimates_11182020.csv', data.table = F) %>%  
  filter(!(Location %in% c('Virgin Islands', 'United States', 'New York City'))) %>%
  filter(as.Date(Sunday, format = '%m/%d/%Y') <= as.Date('4/4/2020',format = '%m/%d/%Y')) %>%
  select(everything(), GLEAM_adj = `Infectious AR median`, -IDEA, - `Infectious AR mean`, -`Infectious AR 2.5% Quantile`, -`Infectious AR 97.5% Quantile`)
 # select(everything(), GLEAM_adj = `Infectious AR median`, -mMAP_adj, -mMAP_IFR065, mMAP_adj = mMAP_adj_IFR065, - `Infectious AR mean`, -`Infectious AR 2.5% Quantile`, -`Infectious AR 97.5% Quantile`)

# create the df of NA's to remove
NA_df = df %>% 
  group_by(Location) %>%
  #select(idea_adj:mMAP) %>% 
  summarize_at(vars(Virology_adj:GLEAM_adj), funs(sum(is.na(.)))) %>%
  tidyr::gather(method, N_NA, Virology_adj:GLEAM_adj) %>%
  filter(N_NA >= 4)

# convert to ratios, convert to long, filter out NAs
df = df %>% 
  group_by(Location) %>%
  summarize_at(vars(Virology_adj:cases), funs(sum(.,na.rm=T))) %>%
  mutate_at(vars(Virology_adj:GLEAM_adj), funs(./cases)) %>%
  select(-cases) %>%
  tidyr::gather(method, N, Virology_adj:GLEAM_adj)  %>%
  na.omit %>%
  anti_join(NA_df, by = c('Location','method'))

# get the min and max across rows (yes there must be a dplyr way to do this)
tmp = df %>% 
  tidyr::spread(method, N) 
tmp$min = apply(tmp[,-1], 1, function(xx) min(xx, na.rm = T))
tmp$max = apply(tmp[,-1], 1, function(xx) max(xx, na.rm = T))
tmp$median = apply(tmp[,-1], 1, function(xx) median(xx, na.rm = T))
df = tidyr::gather(tmp, method, N, covid_scale:median)

# remove IDEA
df = df %>% filter(!(method %in% c('IDEA','IDEA_adj')))

# rename!
df$method = gsub('Virology','div-Vir',df$method)
df$method = gsub('Hist','div-Hist',df$method)
df$method = gsub('covid_scale','COVID Scaling',df$method)
df$method = gsub('_unadj','',df$method)

# box_plot!
df_adj = df[grep('adj',df$method),]
df_un = df[-grep('adj|min|max|median',df$method),]
df_minMax = df[grep('min|max|median',df$method),]

df_adj$method = factor(gsub('_adj','', df_adj$method), 
                       levels = c('div-Hist','div-Vir','COVID Scaling','mMAP','GLEAM'))
df_un$method = factor(gsub('_n','', df_un$method), 
                      levels = c('div-Hist','div-Vir','COVID Scaling','mMAP'))
df_minMax$method = factor(df_minMax$method, levels = c('min','median','max'))

p1 = ggplot(df_adj, aes(x=method, y = N)) + 
  geom_boxplot() + 
  coord_trans(y = 'log1p') + 
  #scale_y_continuous(trans = log1p_trans(), breaks = c(1, 10, 100, 500)) +
  scale_y_continuous(limits = c(0,350), breaks = c(1, 10, 20, 50, 100, 200)) +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        panel.grid.major.y = element_line(colour = 'grey'),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(hjust = 1, vjust = -5)) + 
  labs(y='estimated/reported cases',
       title = 'Adjusted')

p2 = ggplot(df_un, aes(x=method, y = N)) + 
  geom_boxplot() +  
  coord_trans(y = 'log1p') + 
  scale_y_continuous(limits = c(0,350), breaks = c(1, 10, 20, 50, 100, 200)) +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_line(colour = 'grey'),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(hjust = 1, vjust = -5)) + 
  labs(title = 'Unadjusted')

p3 = ggplot(df_minMax, aes(x=method, y = N)) + 
  geom_boxplot() +  
  coord_trans(y = 'log1p') + 
  scale_y_continuous(limits = c(0,350), breaks = c(1, 10, 20, 50, 100, 200)) +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_line(colour = 'grey'),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(hjust = 0.1, vjust = -5)) + 
  labs(title = 'All')

plot_grid(p1, p2, p3, nrow = 1, align = 'h', rel_widths = c(4.8,4,3.2))
ggsave(filename = 'results/boxplot_final_11182020.png', width = 8.5,
       height = 3, unit = 'in')
