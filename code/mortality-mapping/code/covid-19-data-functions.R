# get the moving average of a vector
moving_average <- function(vec, dist = 3){
  v2 = c()
  for(i in 1:length(vec)){
    v2[i] = mean(vec[(max(1, i - dist)):(min(length(vec), i + dist))])
  }
  return(v2)
}

# convert cumulative to incident cases/deaths/hospitalizations
cumulative_to_incident <- function(cases){
  # convert to new cases rather than cumulative sum
  for(i in length(cases):2){cases[i] = cases[i] - cases[i-1]}
  
  # there shouldnt be any negative daily cases
  if(any(cases < 0 )){
    test = cases
    while(any(cases < 0)){
      for(ind in which(cases < 0)){
        # subtract from recent data
        cases[ind - 1] = cases[ind - 1] + cases[ind]
        cases[ind] = 0
      }
    }
  }
  return(cases)
}

# gets the case data for a given state
get_case_data <- function(D, state, country = NULL, abb = NULL, add_starting_zeros = T, smoothing_dist = 0){
  # the columns with data by date
  date_cols = grep('[0-9]/[0-9]',colnames(D))
  
  if(!is.null(country)){
    D = D[D$`Country/Region` == country,]
    cases = colSums(D[,date_cols])
  }else{
    # get the state data
    D_state = D[D$state == state,]
    
    # get state sums
    sum_state = as.numeric(D_state[1, date_cols])
    
    if(!is.null(abb)){
      # get the county data (fixing DC abbr)
      if(abb == 'DC'){
        D_county = D[grepl('D.C.', D$state),] 
      }else{
        D_county = D[grepl(paste0(', ', abb), D$state),]
      }
      
      # get county sums
      sum_county = colSums(D_county[,date_cols])
      
      # join the sums together
      sum_state = sapply(1:length(date_cols), function(ii) max(sum_county[ii], sum_state[ii]))
      
    }
    
    names(sum_state) = colnames(D)[date_cols]
    
    # # skip states with no cases
    # if(sum(sum_state) == 0){
    #   print(sprintf('no cases in %s. Skipping', state))
    # }
    
    # convert to new cases (not cumsum)
    cases = sum_state
  }
  
  # convert to new cases rather than cumulative sum
  for(i in length(cases):2){cases[i] = cases[i] - cases[i-1]}
  
  # there shouldnt be any negative daily cases
  if(any(cases < 0 )){
    test = cases
    while(any(cases < 0)){
      for(ind in which(cases < 0)){
        # subtract from recent data
        cases[ind - 1] = cases[ind - 1] + cases[ind]
        cases[ind] = 0
      }
    }
  }
  
  # convert format of names
  names(cases) = as.Date(names(cases), format = '%m/%d/%y')
  
  # smooth them!
  if(smoothing_dist > 0){
    cases[] = moving_average(cases, 3)
  }
  
  if(!add_starting_zeros){
    cases = cases[cases>0]
  }
  
  return(cases)
}

# gets the case data for a given state
get_county_case_data <- function(D, county, state, add_starting_zeros = T, smoothing_dist = 0){
  # the columns with data by date
  date_cols = grep('[0-9]/[0-9]',colnames(D))
  
  # get the coutny/state data
  D_state = D[D$county == county & D$state == state,]
  if(nrow(D_state) != 1){browser()}
  
  # get state sums
  sum_state = as.numeric(D_state[1, date_cols])
  
  names(sum_state) = colnames(D)[date_cols]

  # convert to new cases (not cumsum)
  cases = sum_state
  
  # convert to new cases rather than cumulative sum
  for(i in length(cases):2){cases[i] = cases[i] - cases[i-1]}
  
  # there shouldnt be any negative daily cases
  if(any(cases < 0 )){
    test = cases
    while(any(cases < 0)){
      for(ind in which(cases < 0)){
        # subtract from recent data
        cases[ind - 1] = cases[ind - 1] + cases[ind]
        cases[ind] = 0
      }
    }
  }
  
  # convert format of names
  names(cases) = as.Date(names(cases), format = '%m/%d/%y')
  
  # smooth them!
  if(smoothing_dist > 0){
    cases[] = moving_average(cases, 3)
  }
  
  if(!add_starting_zeros){
    cases = cases[cases>0]
  }
  
  return(cases)
}

# back-simulates the onset of illness given reported cases
backtrack_cases <- function(cases, testing_delay, expose_rate){
  # create a long list of cases
  long = as.Date(rep(names(cases), cases))
  
  # backtrack the cases
  long = long - (sapply(long, function(xx) min(rgeom(1, expose_rate), ceiling(2/expose_rate))) + 1 + testing_delay)
  
  # create the value to return
  tt = table(long)
  new_cases = as.numeric(tt); names(new_cases) = names(tt)
  
  return(new_cases)
}

# simulate deaths given cases
simulate_deaths <- function(cases, p_death, IFR = .013){
  # initialize it
  deaths = c()
  
  # for each day
  for(i in 1:length(cases)){
    # get current date
    date = as.Date(names(cases)[i])
    
    # simulate number of people who will die
    num_deaths = rbinom(1, cases[i], IFR)
    
    # simulate the dates of these deaths
    deaths = c(deaths, as.character(date + sample(length(p_death), num_deaths, p_death, replace = T)))
  }
  
  # convert the list of deaths to a workable vector
  tt = table(deaths)
  deaths = as.integer(tt); names(deaths) = names(tt)
  min_date = min(as.Date(names(deaths)))
  max_date = max(as.Date(names(deaths)))
  date_seq = seq(min_date, max_date, by = 1)
  missing_dates = setdiff(as.character(date_seq), names(deaths))
  extra_vec = rep(0, length(missing_dates)); names(extra_vec) = missing_dates
  deaths = c(deaths, extra_vec)
  deaths = deaths[order(as.Date(names(deaths)))]
  
  return(deaths)
}

# Creates a transition matrix that multiplies a vector of cases to create a vector of deaths
make_case_to_death_transition <- function(death_dates, case_dates, death_df){
  P = matrix(0, nrow = length(death_dates), ncol = length(case_dates))
  rownames(P) = as.character(death_dates)
  colnames(P) = as.character(case_dates)
  
  # cycle through each row to update
  for(i in 1:nrow(P)){
    # get the case to death differences
    date_diff = as.Date(death_dates[i]) - case_dates
    
    # match these differences with transition vector
    matid = match(date_diff, death_df$days_from_onset)
    
    # update P with the probabilities (will be NA for ones with no transition)
    P[i,] = death_df$p[matid]
  }
  
  # replace NAs with 0
  P[is.na(P)] <- 0
  
  return(P)
}

# run mMAP using matrix operations. Much faster than the original. (this is gonzalo's function for running mMAP)
deaths_to_cases <- function(death, case_to_death, Cvec0 = NULL, max_iter=500, stopping_criteria = 'chi', IFR = 0.013, truncate=TRUE, verbose = T){
  if(!(stopping_criteria %in% c('chi', 'chi-2'))){
    error('havent coded in stopping criteria that is not chi-squared')
  }
  
  # # adding in dates with no deaths
  # dates= seq.Date(min(death$date), max(death$date), by="1 day") %>% as_tibble() %>% rename(date=value)
  # death = dates %>% full_join(death) %>% mutate(death = replace(death, which(is.na(death)),0))
  
  min_date = min(death$date[death$death>0])
  death = death %>% filter(date>=min_date)
  ncols=length(death$date)+length(case_to_death)
  nrows = length(death$date)
  dates_cases = seq.Date(min(death$date)-length(case_to_death),max(death$date),by='1 day')
  
  if(!is.null(Cvec0)){
    Cvec = Cvec0
  }else{
    Cvec = matrix(1, nrow=ncols-1, ncol=1)
  }
  
  deathvec = death$death
  
  #translations of the hospitalization to death distribution
  ma=matrix(0, nrow=nrows,ncol=ncols)
  for(i in seq(1,dim(ma)[1])){
    for(j in seq(i+1,min(i+50,dim(ma)[2]))){
      ma[i,j-1] = case_to_death[50-(j-i)+1]
    }
  }
  
  ma = ma[,1:ncols-1]
  
  right_censor_scale = rev(cumsum(rev(ma[nrows,])))
  right_censor_scale[1:(length(right_censor_scale)-50)]=1
  right_censor_scale[is.na(right_censor_scale)]=1
  
  if(truncate==FALSE){
    right_censor_scale[1:length(right_censor_scale)]=1
  }
  
  # initialize the chi-squared vector
  chi_vec = c()
  
  # predicted deaths initially
  prod1=ma %*% Cvec
  
  for(k in seq(1,max_iter)){
    
    # reported deaths divided by predicted deaths
    norm_vec = t(deathvec/t(prod1))
    
    # sum of transition probabilities and normalizing vec
    prod2 = t(ma) %*% norm_vec
    
    # next iteration of Cvec - updating by transition probabilities, normalizing, and right censor scale
    Cvec=  Cvec*(prod2)/(right_censor_scale)
    
    #browser()
    
    # predicted deaths
    prod1=ma %*% Cvec
    
    # the original (incorrect) way of getting the chi-squared value.
    chi = mean((prod1 - deathvec)^2/prod1)  

    # if chi is NaN (meaning we divided by 0), set to 1000
    #if(is.nan(chi) | is.infinite(chi)){chi = 1000}
    
    chi_vec = c(chi_vec, chi)
    
    # update the minimum chi-squared value and save this Cvec if it is the minimum
    if(k == 1){
      min_chi = chi
      argmin_k = 1
    }else{
      if(chi < min_chi){
        min_chi = chi
        argmin_k = k
        min_chi_vec = Cvec
      }
    }
    
    # if chi-squared is less than 1, stop the iterations
    
    if(chi < 1){
      if(verbose){print(sprintf('chi-squared stat = %s after %s iterations', chi, k))}
      min_chi_vec = Cvec
      break
    }
    if(stopping_criteria == 'chi-2' & k > 1){
      if(abs(chi_vec[k] - chi_vec[k-1])/chi_vec[k-1] < .1){
        #browser()
        if(verbose){print(sprintf('chi-squared stat = %s after %s iterations BC OF <10 percent DECREASE (%s, %s)', chi, k, chi, chi_vec[k-1]))}
        min_chi_vec = Cvec
        break
      }
    }
  }
  if(min_chi > 1){
    if(verbose){print(sprintf('max iterations reached. Minimum chi-squared of %s at iteration %s', min_chi, argmin_k))}
  }
  
  # set the vector estimate to either the vector where chi first < 1 or the minimum chi-squared value
  # and scale it by IFR
  vec_estimate = min_chi_vec[,]/IFR
  death_onsets = tibble(date=dates_cases[1:(length(dates_cases)-1)], daily_cases_estimate=vec_estimate)
  
  return_lst = list(cases = death_onsets, right_censor = right_censor_scale, chi_values = chi_vec)
  return(return_lst)
}

# read in the NYT data
get_time_series_data_NYT <- function(dir = 'data/time_series_data/', type = 'states'){
  files = paste0(dir, dir(dir))
  
  # get the prefix of the file name
  prefix = paste0('NYT_time_series_US_', type)
  
  # only keep the ones with the specified prefix
  files = grep(prefix, files, value = T)
  
  # get the date value 
  dates = sapply(files, function(xx) stringr::str_match(xx, sprintf('%s_(.*?).csv',type))[,2])
  
  # get the most recent date
  ind = which.max(as.Date(dates, format = '%m-%d-%Y'))
  
  # pull data
  D = fread(files[ind], data.table = F)
  
  # convert date format to how it is in the JHU data
  D$date = sapply(as.character(D$date), function(xx) paste(c(strsplit(xx,'-')[[1]][2:3],'20'), collapse='/'))

  if(type == 'states'){
    # get cases 
    cases = D[,c('date','state','cases')]
    
    # cast cases to wide and adjust column names
    cases = reshape(cases, idvar = 'state', v.names = 'cases', timevar = 'date', direction = 'wide')
    cases[is.na(cases)] = 0
    colnames(cases) = gsub('cases.','',colnames(cases))
    
    # make the final US row
    cases[nrow(cases)+1,] = NA
    cases$state[nrow(cases)] = 'United States'
    cases[nrow(cases),-1] = colSums(cases[1:(nrow(cases)-1),-1])
    
    # get deaths
    deaths = D[,c('date','state','deaths')]
    
    # cast deaths to wide and adjust column names
    deaths = reshape(deaths, idvar = 'state', v.names = 'deaths', timevar = 'date', direction = 'wide')
    deaths[is.na(deaths)] = 0
    colnames(deaths) = gsub('deaths.','',colnames(deaths))
    
    # make the final US row
    deaths[nrow(deaths)+1,] = NA
    deaths$state[nrow(deaths)] = 'United States'
    deaths[nrow(deaths),-1] = colSums(deaths[1:(nrow(deaths)-1),-1])
  }else{
    # county-level
    
    # merge county and state names
    D$county_state = paste0(D$county, '_', D$state)
    
    # get cases 
    cases = D[,c('date','county_state','cases')]
    
    # cast cases to wide and adjust column names
    cases = reshape(cases, idvar = 'county_state', v.names = 'cases', timevar = 'date', direction = 'wide')
    cases[is.na(cases)] = 0
    colnames(cases) = gsub('cases.','',colnames(cases))
    
    # get deaths
    deaths = D[,c('date','county_state','deaths')]
    
    # cast deaths to wide and adjust column names
    deaths = reshape(deaths, idvar = 'county_state', v.names = 'deaths', timevar = 'date', direction = 'wide')
    deaths[is.na(deaths)] = 0
    colnames(deaths) = gsub('deaths.','',colnames(deaths))
    
    # split up the states and counties
    cases$state = sapply(cases$county_state, function(xx) strsplit(xx, '_')[[1]][2])
    cases$county = sapply(cases$county_state, function(xx) strsplit(xx, '_')[[1]][1])
    cases$county_state = NULL
    cases = cases[,c(ncol(cases)-1,ncol(cases),1:(ncol(cases)-2))]
    
    deaths$state = sapply(deaths$county_state, function(xx) strsplit(xx, '_')[[1]][2])
    deaths$county = sapply(deaths$county_state, function(xx) strsplit(xx, '_')[[1]][1])
    deaths$county_state = NULL
    deaths = deaths[,c(ncol(deaths)-1,ncol(deaths),1:(ncol(deaths)-2))]
  }

  # return it all!
  return_list = list(cases, deaths, date = dates[ind])
  names(return_list) = c('cases','deaths','date')
  return(return_list)
}

# pulls the most recent time series data and the date of that data
get_time_series_data_JHU <- function(dir = 'data/time_series_data/', type = 'confirmed_US', US_only = T){
  # get the files in the data folder
  files = paste0(dir, dir(dir))
  
  # get the prefix of the file name
  prefix = paste0('JHU_time_series_', type)
  
  # only keep the ones with the specified prefix
  files = grep(prefix, files, value = T)
  
  # get the date value 
  dates = sapply(files, function(xx) stringr::str_match(xx, sprintf('%s_(.*?).csv',type))[,2])
  
  # get the most recent date
  ind = which.max(as.Date(dates, format = '%m-%d-%Y'))
  
  # pull data
  D = fread(files[ind], data.table = F)
  
  if(US_only){
    # only keep US
    D = D[D$`Country/Region` == 'US',]
    
    # create the final US sum column
    D[nrow(D)+1,] = rep(NA, ncol(D))
    D[nrow(D), 1] = 'United States'; D[nrow(D), 2] = 'US'
    D[nrow(D), 5:ncol(D)] = colSums(D[1:(nrow(D)-1),5:ncol(D)])
  }
  
  # get date to return
  date = dates[ind]
  return_list = list(D, date)
  
  # return it!
  return(return_list)
  
}

# return onset-death distribution
return_onset_death_df_hosp <- function(return_all = T, return_beginning = F){
  median = 9.1
  mean = 13
  mu = log(median)
  sigma = sqrt(2*(log(mean) - mu))
  #mu = log(13.2)
  #sd <- sqrt(2*(log(14.5) - mu))
  cum_prob = plnorm(1:50, mu, sigma)
  prob = cum_prob
  for(i in length(prob):2){prob[i] = prob[i] - prob[i-1]}
  df = data.frame(days_from_onset = 1:50, p = prob, p_cum = cum_prob)
  
  if(return_beginning){
    df = df[1:max(which(df$p>.001)),]
    return(df)
  }else if(return_all){
    return(df)
  }else{
    # remove small values
    df = df[df$p > .001,]
    return(df)
  }
}

# return onset-death distribution
return_onset_death_df <- function(return_all = T, return_beginning = F){
  median = 17.1
  mean = 20.2
  mu = log(median)
  sigma = sqrt(2*(log(mean) - mu))
  #mu = log(13.2)
  #sd <- sqrt(2*(log(14.5) - mu))
  cum_prob = plnorm(1:50, mu, sigma)
  prob = cum_prob
  for(i in length(prob):2){prob[i] = prob[i] - prob[i-1]}
  df = data.frame(days_from_onset = 1:50, p = prob, p_cum = cum_prob)
  
  if(return_beginning){
    df = df[1:max(which(df$p>.001)),]
    return(df)
  }else if(return_all){
    return(df)
  }else{
    # remove small values
    df = df[df$p > .001,]
    return(df)
  }
}

# getting the positive and negative tests for each state by day
get_testing_data <- function(dir = 'data/testing_data/'){
  # get the files in the data folder
  files = paste0(dir, dir(dir))
  
  # only keep the ones with the specified prefix
  files = grep('states-daily-testing', files, value = T)
  
  # get the date value 
  dates = sapply(files, function(xx) stringr::str_match(xx, 'testing-(.*?).csv')[,2])
  
  # get the most recent date
  ind = which.max(as.Date(dates, format = '%m-%d-%Y'))
  
  # read in data - oo they also have the hospitalized numbers!
  D_1 = fread(sprintf('data/testing_data/states-daily-testing-%s.csv', dates[ind]))
  D_2 = fread(sprintf('data/testing_data/us-daily-testing-%s.csv', dates[ind]))
  
  # fix US data and merge
  D_2$states = NULL
  D_2$state = 'US'
  
  # merge data
  D_1$fips = NULL
  D_full = rbind(D_1, D_2)
  
  # only keep columns we care about
  D_full = D_full[,.(date, state, positive, negative, hospitalized)]
  
  # fixing date type
  D_full$date = as.Date(as.character(D_full$date), format = '%Y%m%d')
  
  # order data
  D_full = D_full[order(D_full$date),]
  
  # replace NA's with 0
  D_full[is.na(D_full)] = 0
  
  # fix date
  D_daily = NULL
  for(state in unique(D_full$state)){
    # getting the state-specific data set
    ind = which(D_full$state == state)
    tmp = D_full[ind,]
    
    # getting daily values from cumulative
    for(i in nrow(tmp):2){
      tmp$positive[i] = tmp$positive[i] - tmp$positive[i-1]
      tmp$negative[i] = tmp$negative[i] - tmp$negative[i-1]
    }
    # fixing negative values lazily
    if(any(tmp$positive < 0) | any(tmp$negative < 0)){
      print(state)
      # print(tmp)
      tmp$positive[tmp$positive < 0] = 0
      tmp$negative[tmp$negative < 0] = 0
    }
    
    # update daily values
    D_daily = rbind(D_daily, tmp)
  }
  
  # make list to return
  lst = list(cumulative = D_full, daily = D_daily)
  return(lst)
}

# getting global testing data
get_global_testing_data = function(){
  tests_global = fread('data/testing_data/tests-vs-confirmed-cases-global-covid-19.csv', data.table=F)
  tests_global = tests_global[!is.na(tests_global$`Total COVID-19 tests`),]
  
  # 59 is March 20
  tmp_merger = data.frame(Year = 1:70, date = as.Date(NA))
  tmp_merger$date = as.Date('2020-03-20') + 1:70 - 59
  tests_global = merge(tests_global, tmp_merger)
  
  # remove NAs
  tests_global = na.omit(tests_global)
  
  # fix names to match JHU data
  tests_global$Entity[tests_global$Entity == 'Taiwan'] = 'Taiwan*'
  tests_global$Entity[tests_global$Entity == 'Czech Republic'] = 'Czechia'
  tests_global$Entity[tests_global$Entity == 'South Korea'] = 'Korea, South'
  tests_global$Entity[tests_global$Entity == 'United States'] = 'US'
  
  # remove unnecessary columns
  tests_global$Year = tests_global$Code = NULL
  
  # reorder and rename
  tests_global = tests_global[,c(1,4,2,3)]
  colnames(tests_global) = c('country','date','tests','cases')
  
  # return positive and negative tests
  tests_global$positive = tests_global$cases
  tests_global$negative = tests_global$tests - tests_global$cases
  
  # actually just remove cases
  tests_global$cases = NULL
  
  return(tests_global)
}

### using deaths, back-predict the cases. This is the mMAP wrapper
#
# remove_last = number of recent days to remove from prediction
# level = country or 'US'
# min_deaths = minimum number of deaths by covid to required to run mMAP
# asymptomatic_scale_down = what proportion of asymptomatics do we want to exclude (this adjust the fatality rato)
# IFR_scale = how much to scale the infection fatality ratio by
# IFR = What should the IFR be? SHould be NULL for doing state estimates. If this and IFR_scale are not null, this won't work
# add_pneumonia = whether we want to add pneumonia deaths to reported deaths
# value_per_million = whether we want the results in terms of cases per million
# stopping_criteria = 'chi' for chi-squared or 'tol' for tolerance in relative change in iterations
# smoothing_dist = number of days smoothing deaths on both sides (so 3 would be smoothing by week)
# cutoff_date = date to cut off the data at if we want to train using less data.

get_death_predicted_case_data <- function(remove_last = 0, level = 'US', min_deaths = 10, asymptomatic_scale_down = 0.4, IFR_scale = 1,  IFR_constant = NULL, add_pneumonia = F,  value_per_million = F, stopping_criteria = 'chi-2', smoothing_dist = 3, cutoff_date = NULL, verbose = T){
  if(!is.null(IFR_constant) & !(is.null(IFR_scale))){
    stop('please input just IFR or IFR scale')
  }
  if(level == 'country'){
    warning('DEAL WITH HONG KONG LATER')
    D_list = get_time_series_data_JHU(type = 'confirmed_global', US_only = F)
    D_death_list = get_time_series_data_JHU(type = 'deaths_global', US_only = F)
    
    # get the global case and death data
    D = D_list[[1]]; date = D_list[[2]]
    D_death = D_death_list[[1]]; date_death = D_death_list[[2]]
  }else if(level == 'US'){
    D_list = get_time_series_data_NYT(type = 'states')
    D = D_list[[1]]; D_death = D_list[[2]]; date = D_list[[3]]
    
    # getting NYC county
    D_list = get_time_series_data_NYT(type = 'counties')
    D_2 = D_list[[1]]; D_death_2 = D_list[[2]]; date_2 = D_list[[3]]
    
    if(date != date_2){
      stop('county and state dates are different')
    }

    load('data/demographic_data/state_hosp_ICU_IFR_rates.RData')
    state_abb = fread('data/demographic_data/state_abbreviations.csv', data.table = F, select = c(1,3))
    state_pop = fread('data/demographic_data/census-state-population.csv', data.table = F)
  }else{
    stop('please input level = US or country')
  }
  if(add_pneumonia){
    pneum_deaths = fread('data/excess_deaths/pneum_inf_excess_2020.csv', data.table = F)
    #pneum_deaths$weekstart = as.Date(pneum_deaths$weekstart, format = '%m/%d/%Y')
    max_pneum_date = max(pneum_deaths$Week.Ending.Date)
  }
  
  # get the onset to death distribution
  transitions = return_onset_death_df()

  if(smoothing_dist > 0 & verbose){print(sprintf('smoothing death counts with a moving average of %s days on both sides', smoothing_dist))}
  
  if(level == 'country'){
    loc_list = D_death$`Country/Region`[D_death[,ncol(D_death)] >= min_deaths]
  }else{
    loc_list = D_death$state[D_death[,ncol(D_death)] >= min_deaths]
    
    # only keeping state-level data
    loc_list = loc_list[!grepl(', ', loc_list)]
    
    # removing cruise ships, Guam, and the Virgin Islands
    loc_list = loc_list[!grepl('Princess|Guam|(Virgin Islands)|(Northern Mariana Islands)|(American Samoa)',loc_list)]
  }
  
  # initialize
  df = NULL
  chi_list = list()
  death_list = list()
  IFR_list = list()
  
  # back-predict cases by deaths by location
  for(loc in loc_list){
    if(verbose){print(loc)}
    tryCatch({
      if(loc == 'Hong Kong'){
        # get case data
        cases = get_case_data(D, state = loc, smoothing_dist = smoothing_dist)
        deaths = get_case_data(D_death, state = loc, smoothing_dist = smoothing_dist)
      }else if(level == 'country'){
        # get case data
        cases = get_case_data(D, country = loc, smoothing_dist = smoothing_dist)
        deaths = get_case_data(D_death, country = loc, smoothing_dist = smoothing_dist)
      }else{ # (level is state here)
        if(loc == 'Guam'){browser()}
        abb = state_abb$Code[state_abb$State == loc]
        # get case data
        cases = get_case_data(D, state = loc, abb = abb, smoothing_dist = smoothing_dist)
        deaths = get_case_data(D_death, state = loc, abb = abb, smoothing_dist = smoothing_dist)
        
        # add in pneumonia deaths to 
        if(add_pneumonia){
          if(loc %in% pneum_deaths$location){
            backup = deaths
            # pull out the location data frame
            tmp = pneum_deaths[pneum_deaths$location == loc,]
            for(i in 1:nrow(tmp)){
              #weekstart = tmp$weekstart[i]
              weekend = tmp$Week.Ending.Date[i]
              if(weekend < max_pneum_date){
                # find the dates in this week
                # ind = which(as.Date(names(deaths)) >= weekstart &
                #               as.Date(names(deaths)) < weekstart + 7)
                ind = which(as.Date(names(deaths)) > weekend - 7 &
                  as.Date(names(deaths)) <= weekend)
              }else{
                
                # get all the days after the last pneumonia week
                ind = which(as.Date(names(deaths)) > weekend)
              }
              # add in unaccounted pneumonia deaths
              if(!is.na(tmp$excess_deaths[i]) & length(ind > 0)){
                deaths[ind] = deaths[ind] + tmp$excess_deaths[i]/7
              }
            }
          }else{
            print(sprintf('%s not found in pneumonia deaths. Not updating it', loc))
            #browser()
          }
        }
        
        # get the index of the state in the state list
        ind = which(sapply(state_list, function(xx) xx$state) == loc)
        # not present in the data (Puerto Rico), just take the US average
        if(length(ind) == 0){ind = which(sapply(state_list, function(xx) xx$state) == 'United States')}
        
        # get the state infected fatality rate
        IFR = state_list[[ind]]$IFR
      }
      
      if(!is.null(cutoff_date)){
        deaths = deaths[as.Date(names(deaths)) <= cutoff_date]
      }
      
      # reset the IFR based on how much it is scaled up/down
      if(!is.null(IFR_scale)){
        IFR = IFR*IFR_scale
      
      # using a constant IFR
      }else{
        IFR = IFR_constant
      }
      
      # get the death onsets!
      deaths_2 = data.frame(date = as.Date(names(deaths)), death = deaths)
      
      tmp = deaths_to_cases(deaths_2, transitions$p, IFR = IFR, max_iter = 500, stopping_criteria = stopping_criteria, verbose = verbose)
      
      death_onsets = tmp$cases
      chi_list[[loc]] = min(tmp$chi_values)
      death_list[[loc]] = cumsum(deaths)[length(deaths) - remove_last]
      IFR_list[[loc]] = IFR
      
      # create cumulative estimate    
      death_onsets$cum_cases_estimate = cumsum(death_onsets$daily_cases_estimate)
    
      # add cases back in
      cases = cases[as.Date(names(cases)) %in% death_onsets$date]
      cum_cases = cumsum(cases)
      matid = match(as.Date(names(cases)), death_onsets$date)
      death_onsets$daily_cases = death_onsets$cum_cases = NA
      death_onsets$daily_cases[matid] = cases
      death_onsets$cum_cases[matid] = cum_cases
      
      # add deaths in
      deaths = deaths[as.Date(names(deaths)) %in% death_onsets$date]
      cum_deaths = cumsum(deaths)
      matid = match(as.Date(names(deaths)), death_onsets$date)
      death_onsets$daily_deaths = death_onsets$cum_deaths = NA
      death_onsets$daily_deaths[matid] = deaths
      death_onsets$cum_deaths[matid] = cum_deaths
      
      # add in days
      death_onsets$days_firstcase = 0
      ind = which(death_onsets$cum_cases>0)
      death_onsets$days_firstcase[ind] = 1:length(ind)
      
      if(value_per_million == T){
        # scale the values to be per million people
        scale = 1e6/state_pop$N[state_pop$State == loc]
        cols = setdiff(colnames(death_onsets), c('date','days_firstcase'))
        death_onsets[,cols] = death_onsets[,cols]*scale
      }
      # add in country and store it
      death_onsets$abb = ifelse(level == 'US', abb, NA)
      death_onsets$location = loc
      #if(loc == 'Iowa'){browser()}
      df = rbind(df, death_onsets)
    }, error = function(e){
      print('mMAP did not work. Skipping this location')
    })
  }
  
  # remove most recent week (when doing death-projections, this is unstable)
  if(remove_last > 0){
    df[df$date > (max(df$date)-remove_last), c('daily_cases_estimate', 'cum_cases_estimate')] = NA
    if(verbose){print(sprintf('Estimates are as of %s', max(df$date)-remove_last))}
  }
  
  if(asymptomatic_scale_down > 0){
    print(sprintf('scaling estimates to symptomatic cases based on an estimate of %s percent of infections being asymptomatic', 100*asymptomatic_scale_down))
    df$daily_cases_estimate = df$daily_cases_estimate*(1-asymptomatic_scale_down)
    df$cum_cases_estimate = df$cum_cases_estimate*(1-asymptomatic_scale_down)
  }
  lst = list(df = df, chi_list = chi_list, death_list = death_list, IFR_list = IFR_list)
  return(lst)
}

