# Simulation study to investigate ALM and reference prediction
# Function to run each of the analysis types on simulated data
# Howard Thom 21-December-2019

# V2 uses trycatch() to allow recovery from failure in RP or ALM models
# V3 includes random effects models

# Note .fe subscript is used for both random and fixed study effects models - apologies

models.directory<-"models"



require(R2OpenBUGS)

# The following are the NICE TSD2 1c (RE) and 1d (FE) functions
source(paste0(models.directory,"/independent.baselines.model.R"))



# The following uses random effects on baseline to include single arm studies 
# and keeps RCT network otherwise separate from single-arm studies
#source(paste0(models.directory,"/random.effects.baseline.model.single-arms.3.R"))
# The following use fixed effects with (up to) 3 covariates on baseline to include single arm studies and keeps RCT network otherwise separate from single-arm studies
source(paste0(models.directory,"/random.effects.baseline.model.single-arms.1cov.2.R"))


# The following uses random effects on baseline to include disconnected RCTs 
# but keeps RCT network (with reference) otherwise separate
#source(paste0(models.directory,"/random.effects.baseline.model.disconnected.5.R"))
# The following uses fixed effects with (up to) 3 covariates on baseline 
# to include disconnected RCTs but keeps RCT network (with reference) otherwise separate
source(paste0(models.directory,"/random.effects.baseline.model.disconnected.1.cov.3.R"))

# Random effects random effect on baseline models
model.file.random.effects.baseline.single.arms.1cov <- paste0(models.directory,"/random.effects.baseline.model.single-arms.1cov.re.3.txt")
model.file.random.effects.baseline.disconnected.1cov <- paste0(models.directory,"/random.effects.baseline.model.disconnected.1cov.re.4.txt")


# Plugin estimator models
source(paste0(models.directory,"/plugin.model.single.3.R"))
source(paste0(models.directory,"/plugin.model.disconnected.3.R"))

# Random effects plugin estimator models
model.file.plugin.single.re <- paste0(models.directory,"/plugin.model.single.re.3.txt")
model.file.plugin.disconnected.re <- paste0(models.directory,"/plugin.model.disconnected.re.3.txt")

# Default is that it does fixed effects
analyse.simulation.data<-function(b.data, n.chains = 2, 
	num.sims = 3000 * 2, burn.in= 3000 * 2,
	do.debug = FALSE, random.effects = FALSE)
{
	# Set up a summary matrix to use if the BUGS code fails
	# This ensures a row of NA will be exported for each failure in the results tables
	error.matrix <- matrix(NA, nrow=4, ncol=3, dimnames = list(rep("d",4),c("mean","2.5%","97.5%")))


	# Initial values are shared across all models
	# Some values are simply unused by some models
	inits1<-list(d = c(NA, rep(0.5,b.data$nt-1)), mu=rep(0.5,b.data$ns), sd = 1,
			mu.single = rep(-0.5, b.data$ns.single), m = 0.1, sd.m = 1, beta = 0.1,
			mu.base = rep(-0.5, b.data$ns.base), beta.base = 0.1,
			mu.disc = matrix(-0.5, nrow=b.data$ns.disc, ncol=max(b.data$na.disc)),
			sd = 1) # mu.disc is a vector in ALM and matrix in RP
	inits2<-list(d = c(NA, rep(-0.5,b.data$nt-1)), mu=rep(-0.5,b.data$ns), sd = 0.5,
			mu.single = rep(0.5, b.data$ns.single), m = 0.5, sd.m = 0.5, beta = 0.25,
			mu.base = rep(0.5, b.data$ns.base), beta.base = 0.25,
			mu.disc = matrix(0.5, nrow=b.data$ns.disc, ncol=max(b.data$na.disc)),
			sd = 0.5) # mu.disc is a vector in ALM and matrix in RP
	bugs.inits<-list(inits1,inits2)




	# Standard method (independent baselines)
	if(do.debug)print("RCT only")
	if(!random.effects) {
	  bugs.object.rct.only.fe<- tryCatch(
	    bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("mu","d"),model=model.independent.baseline.fe,
	                                clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=do.debug),
	    error = function(e) {list("summary" = error.matrix)})
	} else {
	  bugs.object.rct.only.fe<- tryCatch(
	    bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("mu","d","sd"),model=model.independent.baseline.re,
	                                clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=2,debug=do.debug),
	    error = function(e) {list("summary" = error.matrix)})
	}
	# If the error matrix has been returned
	if(sum(rowSums(is.na(bugs.object.rct.only.fe$summary))) == 12) {
	  # Set all bugs objects to error matrix
	  # For error correcting need all objects to be NULL
	  bugs.object.re.base.disc.fe <- bugs.object.re.base.single.fe <- 
	    bugs.object.plugin.disc.fe <- bugs.object.plugin.single.fe <- NULL
	  
	  bugs.object.re.base.disc.fe$summary <- bugs.object.re.base.single.fe$summary <-
	    bugs.object.plugin.disc.fe$summary <- bugs.object.plugin.single.fe$summary <- error.matrix
	} else {
	  # Do the usual code
	  # Vectors of RCTs matched to the disconnected RCTs and single-arm studies
	  # Used for plug-in models
	  matched.rct.disc <- rep(NA,b.data$ns.disc)
	  matched.rct.single <- rep(NA, b.data$ns.single)
	  
	  # For ALM, need to extract the plug in for the closest matching study
	  # Need separate data structures for single-arm and disconnected analyses for the plug-in models
	  b.data.disc <- b.data.single <- b.data
	  for(i.disc in 1:b.data$ns.disc)
	  {
	    # Use euclidean distance to choose closest arm match
	    distance <- rep(NA, b.data$ns)
	    for(i.rct in 1:b.data$ns)
	    {
	      # Distance between unweighted average of arms of each trial
	      distance[i.rct] <- dist(rbind(
	        mean(b.data$x.disc[i.disc,1:b.data$na.disc[i.disc]]),
	        mean(b.data$x[i.rct,1:b.data$na[i.rct]])
	      ))
	    }
	    min.rct<-which.min(distance)
	    matched.rct.disc[i.disc]<-min.rct
	  }
	  for(i.single in 1:b.data$ns.single)
	  {
	    # Use euclidean distance to choose closest arm match
	    distance <- rep(NA, b.data$ns)
	    for(i.rct in 1:b.data$ns)
	    {
	      # Distance between unweighted average of arms of each trial
	      distance[i.rct] <- dist(rbind(
	        b.data$x.single[i.single],
	        mean(b.data$x[i.rct,1:b.data$na[i.rct]])
	      ))
	    }
	    min.rct<-which.min(distance)
	    matched.rct.single[i.single]<-min.rct
	  }
	  b.data.single$matched.rct <- matched.rct.single
	  b.data.disc$matched.rct <- matched.rct.disc
	  
	  # Take the mu from the matched RCT
	  # Adding NA at end to avoid confusion between vectors and scalars
	  b.data.single$mu.plugin.mean<-
	    c(bugs.object.rct.only.fe$summary[matched.rct.single,"mean"],NA)
	  b.data.single$mu.plugin.prec<-
	    1/c(bugs.object.rct.only.fe$summary[matched.rct.single,"sd"],NA)^2
	  b.data.disc$mu.plugin.mean<-
	    c(bugs.object.rct.only.fe$summary[matched.rct.disc,"mean"],NA)
	  b.data.disc$mu.plugin.prec<-
	    1/c(bugs.object.rct.only.fe$summary[matched.rct.disc,"sd"],NA)^2
	  
	  
	  # And informative priors for sd.disc from the connected RCTs
	  # This is for both plug-in and reference prediction models
	  if(random.effects) {
	    b.data.disc$sd.connected.mean<-b.data.single$sd.connected.mean<-b.data$sd.connected.mean<-
	      bugs.object.rct.only.fe$summary["sd","mean"]
	    b.data.disc$sd.connected.tau<-b.data.single$sd.connected.tau<-b.data$sd.connected.tau<-
	      1/bugs.object.rct.only.fe$summary["sd","sd"]^2
	  }
	  
	  
	  
	  # For error correcting need all objects to be NULL
	  bugs.object.re.base.disc.fe <- bugs.object.re.base.single.fe <- 
	    bugs.object.plugin.disc.fe <- bugs.object.plugin.single.fe <- NULL
	  
	  # Reference prediction with 1 covariate for disconnected studies
	  if(do.debug)print("RP disc")
	  if(!random.effects) {
	    bugs.object.re.base.disc.fe <- tryCatch(
	      bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("d"),model=model.random.effects.baseline.disc.1cov.fe,
	           clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=do.debug),
	      error = function(e) {list("summary" = error.matrix)})
	  } else {
	    bugs.object.re.base.disc.fe <- tryCatch(
	      bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("d"),model.file= model.file.random.effects.baseline.disconnected.1cov,
	           clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=2,debug=do.debug),
	      error = function(e) {list("summary" = error.matrix)})
	  }
	  
	  # Reference prediction with 1 covariate for single-arm studies
	  if(do.debug)print("RP single")
	  if(!random.effects) {
	    bugs.object.re.base.single.fe <- tryCatch(
	      bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("d"),model=model.random.effects.baseline.single.1cov.fe,
	           clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=2,n.thin=1,debug=do.debug),
	      error = function(e) {list("summary" = error.matrix)})
	  } else {
	    bugs.object.re.base.single.fe <- tryCatch(
	      bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("d"),model.file=model.file.random.effects.baseline.single.arms.1cov,
	           clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=2,n.thin=2,debug=do.debug),
	      error = function(e) {list("summary" = error.matrix)})
	  }
	  
	  
	  # Plug in model disconnected
	  if(do.debug)print("ALM disc")
	  if(!random.effects) {
	    bugs.object.plugin.disc.fe <- tryCatch(
	      bugs(data=b.data.disc,inits=bugs.inits,parameters.to.save=c("d"),model=model.plugin.disc.fe,
	           clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=do.debug),
	      error = function(e) {list("summary" = error.matrix)})
	  } else {
	    bugs.object.plugin.disc.fe <- tryCatch(
	      bugs(data=b.data.disc,inits=bugs.inits,parameters.to.save=c("d"),model.file=model.file.plugin.disconnected.re,
	           clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=2,debug=do.debug),
	      error = function(e) {list("summary" = error.matrix)})
	  }
	  
	  # Plug in model single-arm
	  if(do.debug)print("ALM single")
	  if(!random.effects) {
	    bugs.object.plugin.single.fe <- tryCatch(
	      bugs(data=b.data.single,inits=bugs.inits,parameters.to.save=c("d"),model=model.plugin.single.fe,
	           clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=do.debug),
	      error = function(e) {list("summary" = error.matrix)})  
	  } else {
	    bugs.object.plugin.single.fe <- tryCatch(
	      bugs(data=b.data.single,inits=bugs.inits,parameters.to.save=c("d"),model.file=model.file.plugin.single.re,
	           clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=2,debug=do.debug),
	      error = function(e) {list("summary" = error.matrix)})  
	  }
	  
	  
	} # End if the base model fails
	


	# For bias only need the point estimates (mean)
	# For the coverage need the upper and lower 95% CrI limits

	mean.table <-
		rbind(bugs.object.rct.only.fe$summary[grep("d",rownames(bugs.object.rct.only.fe$summary))[1:4],"mean"],
		bugs.object.plugin.single.fe$summary[grep("d",rownames(bugs.object.plugin.single.fe$summary))[1:4],"mean"],
		bugs.object.plugin.disc.fe$summary[grep("d",rownames(bugs.object.plugin.disc.fe$summary))[1:4],"mean"],
		bugs.object.re.base.single.fe$summary[grep("d",rownames(bugs.object.re.base.single.fe$summary))[1:4],"mean"],
		bugs.object.re.base.disc.fe$summary[grep("d",rownames(bugs.object.re.base.disc.fe$summary))[1:4],"mean"])
	ll.table <-
		rbind(bugs.object.rct.only.fe$summary[grep("d",rownames(bugs.object.rct.only.fe$summary))[1:4],"2.5%"],
		bugs.object.plugin.single.fe$summary[grep("d",rownames(bugs.object.plugin.single.fe$summary))[1:4],"2.5%"],
		bugs.object.plugin.disc.fe$summary[grep("d",rownames(bugs.object.plugin.disc.fe$summary))[1:4],"2.5%"],
		bugs.object.re.base.single.fe$summary[grep("d",rownames(bugs.object.re.base.single.fe$summary))[1:4],"2.5%"],
		bugs.object.re.base.disc.fe$summary[grep("d",rownames(bugs.object.re.base.disc.fe$summary))[1:4],"2.5%"])
	ul.table <-
		rbind(bugs.object.rct.only.fe$summary[grep("d",rownames(bugs.object.rct.only.fe$summary))[1:4],"97.5%"],
		bugs.object.plugin.single.fe$summary[grep("d",rownames(bugs.object.plugin.single.fe$summary))[1:4],"97.5%"],
		bugs.object.plugin.disc.fe$summary[grep("d",rownames(bugs.object.plugin.disc.fe$summary))[1:4],"97.5%"],
		bugs.object.re.base.single.fe$summary[grep("d",rownames(bugs.object.re.base.single.fe$summary))[1:4],"97.5%"],
		bugs.object.re.base.disc.fe$summary[grep("d",rownames(bugs.object.re.base.disc.fe$summary))[1:4],"97.5%"])

	rownames(mean.table) <- rownames(ll.table) <- rownames(ul.table) <- c("RCT Only","ALM single","ALM disc","RP single","RP disc")

	return(list(mean.table = mean.table,
			ll.table = ll.table,
			ul.table = ul.table))
}
