# Adaptation of independent baselines models to use plugin estimators for mu
# Howard Thom 20-August-2019
# The closest matching RCT must be matched externally and the 'mu' taken from an independent baselines model

# V2 uses separate standard deviation for random effecs models in connected and disconnected portions
# V3 removes the random effects in order to use the dnorm(,)I(,) syntax


# Extra data are mu.plugin[] and matched.rct[]. The latter is an indicator for the matched RCT
# Data are same as TSD format: ns, nt, na, r, n, t
# The data on disconnected networks are (as in standard TSD):
# ns.disc, nt.disc, na.disc, r.disc, n.disc, t.disc
# x.disc is matrix of covariates for disconnected RCTs.

# sd.connected.mean and sd.connected.tau are the mean and precision of the sd in the connected components.
# These are used as informative priors on sd.disc in the random effects models



# The independent effects on baseline and fixed treatment effect model
model.plugin.disc.fe<-function()
{
	for(i in 1:ns){ # LOOP THROUGH STUDIES 
		mu[i] ~ dnorm(0,.01) # vague priors for all trial baselines 
		for (k in 1:na[i]) { # LOOP THROUGH ARMS 
			r[i,k] ~ dbin(p[i,k],n[i,k]) # Binomial likelihood 
			logit(p[i,k]) <- mu[i] + delta[i,k]
			delta[i,k]<-d[t[i,k]] - d[t[i,1]] # model for linear predictor 
			rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
			dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) 
			+ (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k]))) #Deviance contribution 
		} 
	resdev.rct[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial 
	} 


	for(i in 1:ns.disc){ # LOOP THROUGH STUDIES 
		mu.plugin[i]~dnorm(mu.plugin.mean[i],mu.plugin.prec[i])
		for (k in 1:na.disc[i]) { # LOOP THROUGH ARMS 
			r.disc[i,k] ~ dbin(p.disc[i,k],n.disc[i,k]) # Binomial likelihood 
			logit(p.disc[i,k]) <- mu.plugin[i] + delta.disc[i,k]
			delta.disc[i,k]<-d[t.disc[i,k]] - d[t[matched.rct[i],1]] # model for linear predictor 
			rhat.disc[i,k] <- p.disc[i,k] * n.disc[i,k] # expected value of the numerators 
			dev.disc[i,k] <- 2 * (r.disc[i,k] * (log(r.disc[i,k])-log(rhat.disc[i,k])) 
			+ (n.disc[i,k]-r.disc[i,k]) * (log(n.disc[i,k]-r.disc[i,k]) - log(n.disc[i,k]-rhat.disc[i,k]))) #Deviance contribution 
		} 
	resdev.disc[i] <- sum(dev.disc[i,1:na.disc[i]]) # summed residual deviance contribution for this trial 
	} 
	

	totresdev.rct<-sum(resdev.rct[])
	totresdev.disc<-sum(resdev.disc[])

	totresdev <- totresdev.rct + totresdev.disc  #Total Residual Deviance 
	d[1]<-0 # treatment effect is zero for reference treatment 
	for (k in 2:nt){ d[k] ~ dnorm(0,.01) } # vague priors for treatment effects 
}