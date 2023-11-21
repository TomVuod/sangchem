#' @export
#'
#' @import lmerTest
linear_model<-function(data, MS_data,dependent_var="mass", size_correction=TRUE,
                       peak_subset=NULL,
                       fix_effs=list(), rand_effs=list(), transformation=identity,
                       species="F. sanguinea",...){
  data=data[data$species==species,]
  #colony after queen removal disregarded
  if (any(is.na(data[,dependent_var]))){
    data<-data[!is.na(data[,dependent_var]),]
    warning("Data included NA values")
  }
  if (length(rand_effs)==0) model_call='lm'
  else model_call='lmer'
  if (("mass" %in% colnames(data))&size_correction)
    data$mass<-data$mass/(data$head_width/10^3)^2
  if (any(is.na(data[,"mass"]))) stop("Size correction generated NA values")
  if (!is.null(peak_subset)) data$mass<-
    data$mass*mapply(function(x) subsample_proportion(peak_subset,x, MS_data),
                     data$chromatogram_ID)
  if (any(is.na(data[,"mass"]))) stop("Peak subsetting generated NA values")
  data[,dependent_var]<-transformation(data[,dependent_var])
  if (any(!is.finite(data[,dependent_var]))) return(NULL)
  if (any(is.na(data[,dependent_var]))) stop("Transformation generated NA values")
  model_formula=paste0(dependent_var,"~",paste(fix_effs,collapse='+'))
  if (model_call=="lmer") {random_formula<-paste(paste("(",rand_effs,")",sep=''),collapse="+")
  model_formula=paste(model_formula,random_formula,sep='+')
  }
  #print(model_formula)
  mod.res<-eval(call(model_call,formula=eval(parse(text=model_formula)),data=data))
  mod.diag<-model_diagnostics(mod.res)
  #truncace excess of data in the default result of summary method
  mod.sum<-summary(mod.res)
  mod.sum$call<-call(model_call,formula=eval(parse(text=model_formula)))
  list(mod.res,mod.sum,mod.diag)
}

#' @import DHARMa
model_diagnostics<-function(model){
  counter=1
  results=c()
  if (class(model)=="lmerModLmerTest"){
    for (test in c(testUniformity,testOutliers,testDispersion,testQuantiles)){
      results[counter]=test(model,plot=FALSE)$p.value
      counter=counter+1
    }
  }
  results[counter]<-shapiro.test(residuals(model))$p.value
  if (length(results)==5)
    names(results)<-c("testUniformity","testOutliers","testDispersion","testQuantiles","ShapiroTest")
  else
    names(results)<-c("ShapiroTest")
  results
}
