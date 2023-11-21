#linear_model(development_data, fix_effs = list("sang_prop"), rand_effs = list())

library(dplyr)

sample_averaged<-function(dev_data, MS_data, CHC_amounts){
  if (nrow(filter(dev_data,species=="F. sanguinea"))==0 | nrow(filter(dev_data,species=="F. fusca"))==0)
    return(NULL)
  average_fusca<-averaged_profile(filter(dev_data,species=="F. fusca")$chromatogram_ID,
                                  MS_data = MS_data,
                                  CHC_amount_data = CHC_amounts,
                                  proportions = FALSE)
  average_sanguinea<-averaged_profile(filter(dev_data,species=="F. sanguinea")$chromatogram_ID,
                                      MS_data = MS_data,
                                      CHC_amount_data = CHC_amounts,
                                      proportions = FALSE)
  if(is.null(average_fusca)|is.null(average_sanguinea)) return(NULL)
  joined<-full_join(average_fusca,average_sanguinea,by="peak_ID")
  joined$abundance.x[is.na(joined$abundance.x)]<-0
  joined$abundance.y[is.na(joined$abundance.y)]<-0
  joined
}

extract_model_parms<-function(data){
  model_res<-lm(abundance.x~abundance.y,data=data)
  #summary(model_lm)$adj.r.squared
  rsq=summary(model_res)$r.squared
  n=nrow(data)
  sd.x=sd(data$abundance.x)
  sang_mean<-mean(data$abundance.y)
  fusca_mean<-mean(data$abundance.x)
  correlation=cor(data$abundance.x,data$abundance.y)
  a=model_res$coeff[2]
  p_val <- ifelse(nrow(model_res$model)>1,summary(model_res)$coeff[2,4], NA)
  data.frame(r.squared=rsq,n=n,sd.x=sd.x,
             sang_mean=sang_mean,fusca_mean=fusca_mean,
             correlation=correlation,a=a, p_val = p_val)
}

transferability<-function(data){
  data<-split(data, data$peak_ID)
  data<-lapply(data, extract_model_parms)
  data<-do.call(rbind, data)
  data$peak_ID<-as.numeric(rownames(data))
  data
}

averaged_profile<-function(chr_IDs,MS_data,
                           CHC_amount_data, proportions=TRUE){
  if(length(chr_IDs)==0) stop("No chromatogram ID values in averaged_profile function")
  if(any(is.na(chr_IDs))) warning("NA chromatogram ID values in averaged_profile function")
  #if (length(chr_IDs)!=2) stop("Two chromatogram ID values should be passed in")
  output<-data.frame()
  for (i in 1:length(chr_IDs)){
    new_data<-relative_amounts(chr_IDs[i], chromatogram_data= MS_data)
    if (is.null(new_data)) next
    if (!proportions){
      mass=filter(CHC_amount_data,chromatogram_ID==chr_IDs[i])$corrected_mass
      if(is.na(mass)) return(NULL)
      new_data$relative_abundance=mass*new_data$relative_abundance
    }
    new_data$chromatogram_ID<-chr_IDs[i]
    output<-rbind(output,new_data)
  }
  if (nrow(output)==0) return(NULL)
  output<-aggregate(output$relative_abundance, list(output$peak_ID,
                                                    output$chromatogram_ID), sum, drop=FALSE)
  output$x[is.na(output$x)]<-0
  output<-rename(output,relative_abundance=x,peak_ID=Group.1,chromatogram_ID=Group.2)
  output<-aggregate(output$relative_abundance, list(output$peak_ID), mean)
  output=rename(output, relative_abundance=x, peak_ID=Group.1)
  if (!proportions){
    output=rename(output, abundance=relative_abundance)
  }
  output
}

relative_amounts<-function(chrID, chromatogram_data, peak_subset=NULL){
  if (!(is.numeric(chrID))) stop("Chromatogram ID should be numeric")
  if (length(chrID)!=1) stop("One chromatogram ID should be passed in")
  if (!(chrID %in% chromatogram_data$chromatogram_ID)) {
    warning(sprintf("Chromatogram ID not found: %s", as.character(chrID)))
    return(NULL)
  }
  filter(chromatogram_data,chromatogram_ID %in% chrID,peak_ID>1) %>%
    select(chromatogram_ID,peak_ID,abundance) -> filtered_data
  if (!is.null(peak_subset)) filtered_data<-filter(filtered_data,peak_ID %in% peak_subset)
  filtered_data<-aggregate(filtered_data$abundance,
                           list(peak_ID=filtered_data$peak_ID),sum)
  filtered_data<-rename(filtered_data,abundance=x)
  filtered_data<-filtered_data[,c("abundance","peak_ID")]
  filtered_data$abundance<-filtered_data$abundance/sum(filtered_data$abundance)
  rename(filtered_data, relative_abundance=abundance)
}

library(dplyr)
data("development_data")
data("CHC_amounts")
data("mass_spectra_data")
development_data%>%filter(callow==0,caste=="worker", chromatogram_ID %in% mass_spectra_data$chromatogram_ID) %>%
  group_by(colony, census_date) %>%
  group_map(~sample_averaged(.x, MS_data = mass_spectra_data, CHC_amounts = CHC_amounts)) -> exchange_data

exchange_data<-do.call(rbind,exchange_data)

# remove peaks with zero values for both species since otherwise it would evelate correlation
exchange_data<-exchange_data[apply(exchange_data,1,function(x) all(x!=0)),]

peak_trans<-transferability(exchange_data)

peak_trans$n_alkane<-peak_trans$peak_ID %in% n_alkanes_ID

peak_trans %>% filter(n>21,correlation <0.2) %>% select(peak_ID) %>%
  unlist()->low_trans_peaks
peak_trans %>% filter(n>21, p_val <= 0.05, correlation > 0.5) %>% select(peak_ID) %>%
  unlist() ->high_trans_peaks

peak_trans %>% filter(n>11,correlation<0.1) %>% select(peak_ID) %>%
  unlist() %>% intersect(discr_sang)->low_trans_peaks_sang
peak_trans %>% filter(n>11,correlation>0.4) %>% select(peak_ID) %>%
  unlist() %>% intersect(discr_sang)->high_trans_peaks_sang

c(low_trans_model_sang)

caclulate_slope_quotients<-function(model_numbers){
  # if differences in compound amounts reflect only the difference in the colony composition
  # then changes in low and high transferability compounds should be proportional between species
  # with the same proportion coefficient
  # Moreover if colony composition does not affect CHC composition but only the worker age
  # the we should expect similar percentage
  res=list(quotients=list(),slopes=list())
  for (i in 1:length(model_numbers)){
    back_transformation=bm_colony_development[[model_numbers[i]]]$arguments$transformation
    if (identical(back_transformation,log))
      back_transformation=exp
    if (identical(back_transformation,sqrt))
      back_transformation=function(x) x^2
    b=bm_colony_development[[model_numbers[i]]]$model_summary$coeff[1]
    a=bm_colony_development[[model_numbers[i]]]$model_summary$coeff[2]
    res$slopes[[i]]=back_transformation(a+b)-back_transformation(b)
  }
  names(res$slopes)=c("low transferability sanguinea","high transferability sanguinea",
                      "low transferability fusca","high transferability fusca")
  res$quotients=list(low=res$slopes[[1]]/res$slopes[[3]],high=res$slopes[[2]]/res$slopes[[4]])
  res
}

caclulate_slope_quotients(c(1,3,6,8))



simulate_extremes<-function(model_numbers){
  res=list()
  for (i in 1:length(model_numbers)){
    back_transformation=bm_colony_development[[model_numbers[i]]]$arguments$transformation
    if (identical(back_transformation,log))
      back_transformation=exp
    if (identical(back_transformation,sqrt))
      back_transformation=function(x) x^2
    res[[i]]=bootMer(bm_colony_development[[model_numbers[i]]]$model,
                     function(x) predict(x,data.frame(sang_prop=c(0,1)), re.form = NA),
                     nsim=10^3)$t
    res[[i]]=back_transformation(res[[i]])
  }
  res
}

simulate_slope_quotients<-function(data){
  slopes<-list(c(),c())
  for(j in 1:2){
    for(i in 1:10^4){
      rand_n1<-sample(1:nrow(data[[j]]),1)
      rand_n2<-sample(1:nrow(data[[j+2]]),1)
      temp1<-data[[j]][rand_n1,2]-data[[j]][rand_n1,1]
      temp2<-data[[j+2]][rand_n2,2]-data[[j+2]][rand_n2,1]
      slopes[[j]][i]<-temp1/temp2
    }
  }
  slopes
}

temp2<-simulate_slope_quotients(temp)
