####GENERAL FUNCTIONS FOR R####
##AUTHOR: SHOM MAZUMDER
##

#packages
require(reshape)
require(dplyr)
require(ggplot2)
require(RCurl)
require(sandwich)
require(multiwayvcov)

#### GGPLOT THEMES ####

#' Shom's Custom ggplot2 themes
#'
#' Custom ggplot theme for general plots
#'
#' @param None
#'
#' @return None
#'
#' @export
theme_shom <- function(){
  theme_minimal() +
  theme(
    text = element_text(family = "Minion Pro", color = "#22211d",size = 12),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
    panel.grid.major = element_line(color = "#ebebe5", size = 0.5),
    panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "none",
    panel.border = element_blank())
}

#' Shom's Custom ggplot2 themes
#'
#' Custom ggplot theme for choropleths
#'
#' @param None
#'
#' @return None
#'
#' @export
theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank()
    )
}

####DATA PREP####

#' Makes a stata codebook
#'
#' @param x Takes in a vector of labelled values
#'
#' @return None
#'
#' @export
stata_codebook <- function(x) {
  cb <- data.frame(attr(x, "var.labels"))
  rownames(cb) <- names(x)
  cb
}

####AUXILLARY####

#' Checks to see if you have an internet connection
#'
#' @param None
#'
#' @return None
#'
#' @examples
#' got_Internet()
#'
#' @export
got_Internet <- function(){
  if (is.character(getURL("www.google.com"))) {
    out <- TRUE
  } else {
    out <- FALSE
  }
}

####REGRESSION####

#' Checks to see if you have an internet connection
#'
#' @param models.list A list of regression models (lm)
#' @param hc Type of robust standard error. Defaults to HC2.
#'
#' @return se.list
#'
#' @examples
#' get_rse()
#'
#' @export
get_rse <- function(models.list,hc='HC2'){
  ##Purpose: this function takes in a list of regression models and outputs a list of robust SEs
  se.list <- lapply(models.list,FUN = function(x){return(sqrt(diag(vcovHC(x,type = hc))))})
  return(se.list)
}

get_cluster_se <- function(list.models,cluster){
  ##Purpose: takes in a list of regression models and a formula for clustering and outputs list of clustered SEs
  clust <- substitute(cluster)
  vcovCL.list <- lapply(list.models,FUN = function(x){return(cluster.vcov(x,cluster = eval(clust)))})
  seCL.list <- lapply(vcovCL.list,FUN = function(x){return(sqrt(diag(x)))})
  return(seCL.list)
}

ch.row <- function(name, yesno) {
  ##Purpose: adds checkmarks to table
  c(name, ifelse(yesno, "$\\checkmark$", ""))
}

####GIS####


####PLOTTING####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##Coefficient Plots##

myCoefPlot <- function(models = models.list,names = model.names,coef.name = coef.name,se = NULL,xlab="Model"){
  #takes in regression models, the name of the coefficient you want to plot, and an optional list of standard errors
  #and returns a ggplot coefficient plot with 95% confidence intervals

  #initialize vector of regression coefficients
  reg.coefs <- as.vector(unlist(lapply(models,FUN = function(x){return(coef(x)[coef.name])})))

  #initialize vector of regression standard errors
  if(is.null(se)){
    #regular standard errors
    se.vector <- unlist(lapply(models,FUN = function(x){return(sqrt(diag(vcov(x)))[coef.name])}))
  }else{
    #user supplied standard errors
    se.vector <- unlist(lapply(se,FUN = function(x){return(x[coef.name])}))
  }

  #generate vector of 95% confidence intervals
  reg.lower.95 <- reg.coefs-qnorm(p = 0.975)*se.vector
  reg.upper.95 <- reg.coefs+qnorm(p = 0.975)*se.vector

  #generate vector of 90% confidence intervals
  reg.lower.90 <- reg.coefs-qnorm(p = 0.95)*se.vector
  reg.upper.90 <- reg.coefs+qnorm(p = 0.95)*se.vector

  #generate dataframe for plotting in ggplot2
  plot.df <- data.frame(cbind(reg.coefs,reg.lower.95,reg.upper.95,reg.lower.90,reg.upper.90,Model=names),row.names = NULL)
  plot.df$reg.coefs <- as.numeric(as.character(plot.df$reg.coefs))
  plot.df$reg.lower.95 <- as.numeric(as.character(plot.df$reg.lower.95))
  plot.df$reg.upper.95 <- as.numeric(as.character(plot.df$reg.upper.95))
  plot.df$reg.lower.90 <- as.numeric(as.character(plot.df$reg.lower.90))
  plot.df$reg.upper.90 <- as.numeric(as.character(plot.df$reg.upper.90))

  #plot coefficients
  reg.coef.plot <- plot.df %>% ggplot(aes(x=Model,y=reg.coefs,group=Model,colour=Model)) +
    geom_point() +
    geom_hline(yintercept = 0,size=0.5) +
    geom_linerange(aes(ymax = reg.upper.95,ymin=reg.lower.95))+
    geom_linerange(aes(ymax = reg.upper.90,ymin=reg.lower.90),size=1)+
    ylab('Estimated Coefficient')+
    xlab(xlab)+
    theme(text = element_text(size=15),legend.position = "none")
  return(reg.coef.plot)
}


####INSTRUMENTAL VARIABLES###

#get first-stage f-stat
getFStat <- function(ivobject){
  sumiv <- summary(ivobject,diagnostics=T)
  f.stat <- as.numeric(sumiv$diagnostics['Weak instruments','statistic'])
  return(round(f.stat,digits = 2))
}

#iv exclusion sensitivity analysis (Conley et al 2012, ReStat)
#code built off of Jim Bisbee's function (http://static1.squarespace.com/static/559c0c95e4b0c021321c7435/t/56d4a18307eaa0ae80a870eb/1456775555253/handout-2-29.pdf)
#exclusion_sens <- function(ivobject,dat,g){
#  gamma <- seq(-1,1,0.1)
#  newY <- dat$dv - g*dat$instrument
#}
#gamma <- seq(-1,1,.25)
#ExclSens <- function(g) {
#  newY <- dat$logpgp95 - g*dat$logem4
#  coef(ivreg(newY~avexpr+f_brit+f_french,~logem4+f_brit+f_french,cbind(dat,newY)))[2]
#}
#sens.coefs <- sapply(gamma,ExclSens)
#names(sens.coefs)<-gamma
#round(sens.coefs,3)

####DIFF IN DIFF####
#function to create parallel trend plot (NEEDS TO BE FIXED)
parallel.trend <- function(dv,upper,lower,df,treatment){#need to generalize this
  dv <- substitute(dv)
  y.max <- substitute(upper)
  y.min <- substitute(lower)
  y.lab <- NULL
  title <- NULL
  parallel.plot <- ggplot(df,aes(x=year,y=eval(dv),group=treatment,color=factor(treatment)))+
    geom_point()+
    geom_line()+
    geom_vline(xintercept = 2009,size=1.5)+
    geom_linerange(aes(ymax=eval(y.max),ymin=eval(y.min)))
  return(parallel.plot)
}

## Two-variable interaction plots in R
## Anton Strezhnev
## 06/17/2013

## interaction_plot_continuous: Plots the marginal effect for one variable interacted with a continuous moderator variable
## Usage
## Required
# model: linear or generalized linear model object returned by lm() or glm() function
# effect: name of the "effect" variable in the interaction (marginal effect plotted on y-axis) - character string
# moderator: name of the moderating variable in the interaction (plotted on x-axis) - character string
# interaction: name of the interaction variable in the model object - character string
## Optional
# varcov: Variance-Covariance matrix - if default, then taken from the model object using vcov()
# minimum: Smallest value of moderator for which a marginal effect is calculated, if "min" then equal to the minimum value of the moderator in the dataset
# maximum: Largest value of moderator for which a marginal effect is calucated, if "max" then equal to the maximum value of the moderator in the dataset
# num_points: Total number of points for which a marginal effect is calculated - increase to make confidence bounds appear smoother
# conf: Size of confidence interval around coefficient estimates - 0-1, default is .95 (95% confidence)
# mean: Mark the mean mediator value by a vertical red line
# median: Mark the median mediator value by a vertical blue line
# alph: Transparency level of the histogram plot - 0-100, decrease to make the histogram more transparent
# rugplot: Include a rug plot of the mediator values below the figure
# histogram: Include a histogram of mediator values behind the figure - only plotted if minimum="min" and maximum="max"
# title: Title of the plot
# xlabel: Label of the X axis
# ylabel: Label of the Y axis
interaction_plot_continuous <- function(model, effect, moderator, interaction, varcov = NULL, minimum="min", maximum="max", incr="default", num_points = 50, conf=.95, mean=FALSE, median=FALSE, alph=80, rugplot=T, histogram=T, title="Marginal effects plot", xlabel="Value of moderator", ylabel="Estimated marginal coefficient"){

  # Define a function to make colors transparent
  makeTransparent<-function(someColor, alpha=alph){
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
  }

  # Extract Variance Covariance matrix
  if (is.null(varcov)){
    covMat = vcov(model)
  }else{
    covMat = varcov
  }

  # Extract the data frame of the model
  mod_frame = model.frame(model)

  # Get coefficients of variables
  beta_1 = model$coefficients[[effect]]
  beta_3 = model$coefficients[[interaction]]

  # Set range of the moderator variable
  # Minimum
  if (minimum == "min"){
    min_val = min(mod_frame[[moderator]])
  }else{
    min_val = minimum
  }
  # Maximum
  if (maximum == "max"){
    max_val = max(mod_frame[[moderator]])
  }else{
    max_val = maximum
  }

  # Check if minimum smaller than maximum
  if (min_val > max_val){
    stop("Error: Minimum moderator value greater than maximum value.")
  }

  # Determine intervals between values of the moderator
  if (incr == "default"){
    increment = (max_val - min_val)/(num_points - 1)
  }else{
    increment = incr
  }

  # Create list of moderator values at which marginal effect is evaluated
  x_2 <- seq(from=min_val, to=max_val, by=increment)

  # Compute marginal effects
  delta_1 = beta_1 + beta_3*x_2

  # Compute variances
  var_1 = covMat[effect,effect] + (x_2^2)*covMat[interaction, interaction] + 2*x_2*covMat[effect, interaction]

  # Standard errors
  se_1 = sqrt(var_1)

  # Upper and lower confidence bounds
  z_score = qnorm(1 - ((1 - conf)/2))
  upper_bound = delta_1 + z_score*se_1
  lower_bound = delta_1 - z_score*se_1

  # Determine the bounds of the graphing area
  max_y = max(upper_bound)
  min_y = min(lower_bound)

  # Make the histogram color
  hist_col = makeTransparent("DodgerBlue")

  # Initialize plotting window
  plot(x=c(), y=c(), ylim=c(min_y, max_y), xlim=c(min_val, max_val), xlab=xlabel, ylab=ylabel, main=title,bty = "n")

  # Plot estimated effects
  lines(y=delta_1, x=x_2, col = "IndianRed",lwd = 3)
  lines(y=upper_bound, x=x_2, lty=2, col = "IndianRed",lwd = 2)
  lines(y=lower_bound, x=x_2, lty=2, col = "IndianRed",lwd = 2)

  # Add a dashed horizontal line for zero
  abline(h=0, lty=3)

  # Shade confidence interval
  #polygon(c(x_2,rev(x_2)),c(lower_bound,rev(upper_bound)),col = adjustcolor("DodgerBlue",alpha.f=0.5), border = FALSE)


  # Add a vertical line at the mean
  if (mean){
    abline(v = mean(mod_frame[[moderator]]), lty=2, col="IndianRed")
  }

  # Add a vertical line at the median
  if (median){
    abline(v = median(mod_frame[[moderator]]), lty=3, col="blue")
  }

  # Add Rug plot
  if (rugplot){
    rug(mod_frame[[moderator]])
  }
  # Add Histogram (Histogram only plots when minimum and maximum are the min/max of the moderator)
  if (histogram & minimum=="min" & maximum=="max"){
    par(new=T)
    hist(mod_frame[[moderator]], axes=F, xlab="", ylab="",main="", border=hist_col, col=hist_col)
  }
}

### INTERACTION PLOTS ####
## interaction_plot_binary: Plots the marginal effect for one variable interacted with a binary variable
## Usage
## Required
# model: linear or generalized linear model object returned by lm() or glm() function
# effect: name of the "effect" variable in the interaction (marginal effect plotted on y-axis) - character string
# moderator: name of the moderating variable in the interaction (plotted on x-axis) - character string - Variable must be binary (0-1)
# interaction: name of the interaction variable in the model object - character string
## Optional
# varcov: Variance-Covariance matrix - if default, then taken from the model object using vcov()
# conf: Size of confidence interval around coefficient estimates - 0-1, default is .95 (95% confidence)
# title: Title of the plot
# xlabel: Label of the X axis
# ylabel: Label of the Y axis
# factor_labels: Labels for each of the two moderator values - default = "0" and "1"
interaction_plot_binary <- function(model, effect, moderator, interaction, varcov="default", conf=.95, title="Marginal effects plot", xlabel="Value of moderator", ylabel="Estimated marginal coefficient", factor_labels=c(0,1)){

  # Extract Variance Covariance matrix
  if (varcov == "default"){
    covMat = vcov(model)
  }else{
    covMat = varcov
  }

  # Extract the data frame of the model
  mod_frame = model.frame(model)

  # Get coefficients of variables
  beta_1 = model$coefficients[[effect]]
  beta_3 = model$coefficients[[interaction]]

  # Create list of moderator values at which marginal effect is evaluated
  x_2 <- c(0,1)

  # Compute marginal effects
  delta_1 = beta_1 + beta_3*x_2

  # Compute variances
  var_1 = covMat[effect,effect] + (x_2^2)*covMat[interaction, interaction] + 2*x_2*covMat[effect, interaction]

  # Standard errors
  se_1 = sqrt(var_1)

  # Upper and lower confidence bounds
  z_score = qnorm(1 - ((1 - conf)/2))
  upper_bound = delta_1 + z_score*se_1
  lower_bound = delta_1 - z_score*se_1

  # Determine the bounds of the graphing area
  max_y = max(upper_bound)
  min_y = min(lower_bound)

  # Initialize plotting window
  plot(x=c(), y=c(), ylim=c(min_y, max_y), xlim=c(-.5, 1.5), xlab=xlabel, ylab=ylabel, main=title, xaxt="n")

  # Plot points of estimated effects
  points(x=x_2, y=delta_1, pch=16)

  # Plot lines of confidence intervals
  lines(x=c(x_2[1], x_2[1]), y=c(upper_bound[1], lower_bound[1]), lty=1)
  points(x=c(x_2[1], x_2[1]), y=c(upper_bound[1], lower_bound[1]), pch=c(25,24), bg="black")
  lines(x=c(x_2[2], x_2[2]), y=c(upper_bound[2], lower_bound[2]), lty=1)
  points(x=c(x_2[2], x_2[2]), y=c(upper_bound[2], lower_bound[2]), pch=c(25,24), bg="black")

  # Label the axis
  axis(side=1, at=c(0,1), labels=factor_labels)

  # Add a dashed horizontal line for zero
  abline(h=0, lty=3)

}

#### TABLES ####
format_and_save_tab <- function(table,path){
  #formats and saves stargazer tables
  table <- gsub("\\{\\*\\}", "\\{\\\\dagger\\}", table)
  table <- gsub("\\{\\*\\*\\}", "\\{\\*\\}", table)
  table <- gsub("\\{\\*\\*\\*\\}", "\\{\\*\\*\\}", table)
  cat(paste(table, collapse = "\n"), "\n",file = path)
}

#### MY REGRESSION DISCONTINUITY FUNCTIONS ####

#' Shom's Custom ggplot2 themes
#'
#' Custom ggplot theme for choropleths
#'
#' @param data A dataframe
#' @param forcing The forcing/running variable as a string
#' @param outcome The outcome varaible as a string
#' @param cutoff The cutoff. Defaults to 0.
#' @param weights Any weights you want to use as a string. Defaults to NULL.
#' @param ylab The label for the y axis.
#' @param xlab The label for the x axis.
#' @param bw The bandwidth. Defaults to the entire support of the forcing variable.
#' @param se Logical for whether to include the standard errors in the plot.
#'
#' @return a list of the plot (in ggplot) and the corresponding binned data
#'
#' @export
plot_rdd_binned <- function(data,forcing,outcome,cutoff = 0,weights = NULL,
                            ylab = "Outcome",xlab = "Running Variable",
                            bw = NULL,se = T){
  #requirements
  require(lazyeval)
  require(ggplot2)
  require(dplyr)

  #check conditions
  if(is.null(data)){
    stop("Please provide a dataframe.")
  }
  if(is.null(forcing)){
    stop("Please provide a forcing variable.")
  }
  if(is.null(outcome)){
    stop("Please provide an outcome variable.")
  }
  if(!is_character(forcing)){
    stop("Forcing variable must be provided as a string.")
  }
  if(!is_character(outcome)){
    stop("Outcome variable must be provided as a string.")
  }
  if(!is.null(weights) & !is_character(weights)){
    stop("Outcome variable must be provided as a string.")
  }
  #check if running variable is discrete (TO-DO)
  #ASSUMING DATA IS DISCRETE FOR NOW
  length.unique.forcing <- nrow(unique(data[,forcing]))
  length.forcing <- nrow(data[,forcing])
  if(is.null(weights)){ #no weights provided
    binned_df <- data %>%
      group_by_(forcing) %>%
      dplyr::summarize_(meanoutcome = interp(~mean(var,na.rm = T),var = as.name(outcome)),
                        size.bin = "n()")

    #interp(~n_distinct(v), v=as.name(uniq.var))
  }else{
    binned_df <- data %>%
      group_by_(forcing) %>%
      dplyr::summarize_(meanoutcome = interp(~weighted.mean(var,w,na.rm = T),
                                             var = as.name(outcome),w = as.name(weights)),
                        size.bin = "n()")
  }

  #set bandwidth and filter data
  if(is.null(bw)){
    bw <- max(data[,forcing]) - min(data[,forcing])
  }
  binned_df <- binned_df %>%
    dplyr::mutate_(absforcing = interp(~abs(var),var = as.name(forcing)),
                   running = interp(~var,var = as.name(forcing))) %>%
    dplyr::filter(absforcing <= bw)

  #add outcome labels to data
  binned_df$outcomelabel <- ylab

  #now plot
  p <- binned_df %>%
    ggplot(aes(x = running, y = meanoutcome,size = size.bin,weight = size.bin)) +
    geom_point(col = "IndianRed3",fill = "grey30",shape = 21) +
    geom_smooth(data = subset(binned_df, running >= cutoff),se = se,color = "dodgerblue") +
    geom_smooth(data = subset(binned_df, running < cutoff),se = se,color = "dodgerblue") +
    geom_vline(xintercept = cutoff) +
    theme_bw() +
    theme(legend.position="none") +
    xlab(xlab) +
    ylab(ylab)

  return(list(
    plot_out = p,
    binned_df = binned_df
  ))

}
