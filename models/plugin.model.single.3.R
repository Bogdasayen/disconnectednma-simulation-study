# Plugin model for single arm studies
# Howard Thom 20-August-2019

# V2 uses separate standard deviation for random effecs models in connected and single-arm portions
# V3 removes the random effects in order to use the dnorm(,)I(,) syntax


# Data are same as TSD format: ns, nt, na, r, n, t
# The data on single arm studies are:
# ns.single : number of single arm studies
# r.single : number of events in single arm studies
# n.single : number of patients in single arm studies
# t.single : treatment in single arm study
# mu.plugin[] is the plugin estimator from RCT with index matched.rct[]

# sd.connected.mean and sd.connected.tau are the mean and precision of the sd in the connected components.
# These are used as informative priors on sd.disc in the random effects models

# The random effects on baseline and fixed treatment effect model, including single arm studies
model.plugin.single.fe<-function()
{
	# Model for RCTs ###################################
	for(i in 1:ns){ # LOOP THROUGH STUDIES 
		mu[i] ~ dnorm(0,0.01) # random effect on baselines
		for (k in 1:na[i]) { # LOOP THROUGH ARMS 
			r[i,k] ~ dbin(p[i,k],n[i,k]) # Binomial likelihood 
			logit(p[i,k]) <- mu[i] + delta[i,k]
			delta[i,k]<-d[t[i,k]] - d[t[i,1]] # model for linear predictor 
			rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
			dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) 
			+ (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k]))) #Deviance contribution 
		} 
	resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial 
	} 

	# Model for single-arm studies ###################################
	for(i in 1:ns.single)
	{
		mu.plugin[i]~dnorm(mu.plugin.mean[i],mu.plugin.prec[i])
		r.single[i]~dbin(p.single[i],n.single[i])
		logit(p.single[i])<-mu.plugin[i]+delta.single[i]
		delta.single[i]<-d[t.single[i]] - d[t[matched.rct[i],1]] # Treatment effect relative to reference

		rhat.single[i] <- p.single[i] * n.single[i] # expected value of the numerators 
		dev.single[i] <- 2 * (r.single[i] * (log(r.single[i])-log(rhat.single[i])) 
		+ (n.single[i]-r.single[i]) * (log(n.single[i]-r.single[i]) - log(n.single[i]-rhat.single[i]))) #Deviance contribution 
	}

	# Calculate deviance ###################################
	totresdev.single<-sum(dev.single[]) # Total residual deviance for single-arm studies
	totresdev.rct <- sum(resdev[]) # Total residual deviance for RCTs
	totresdev<-totresdev.single+totresdev.rct #Total Residual Deviance 
	
	# Priors for remaining parameters
	d[1]<-0 # treatment effect is zero for reference treatment 
	for (k in 2:nt){ d[k] ~ dnorm(0,.01) } # vague priors for treatment effects 
}



