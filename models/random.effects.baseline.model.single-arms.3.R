# Random effects on baseline models with single-arm studies
# Howard Thom 6-October-2017

# V2 uses separate standard deviation for random effecs models in connected and single-arm portions
# V3 removes the random effects in order to use the dnorm(,)I(,) syntax

# Changes from previous versions
# Uses only reference treatment arms to fit the baseline model
# Uses 'cut()' to prevent feedback from single-arm studies to baseline model
# Naively pools single arm and RCT evidence on treatment effects (not sure enough evidence to fit hierarchical model)

# Data are same as TSD format: ns, nt, na, r, n, t
# The data on baseline arms (from RCTs) are
# ns.base : number of baseline arms
# r.base : number of events in baseline arms
# n.base : number of patients in baseline arms
# The data on single arm studies are:
# ns.single : number of single arm studies
# r.single : number of events in single arm studies
# n.single : number of patients in single arm studies
# t.single : treatment in single arm study

# sd.connected.mean and sd.connected.tau are the mean and precision of the sd in the connected components.
# These are used as informative priors on sd.disc in the random effects models

# The random effects on baseline and fixed treatment effect model, including single arm studies
model.random.effects.baseline.single.fe<-function()
{
	# Model for RCTs ###################################
	for(i in 1:ns){ # LOOP THROUGH STUDIES 
		mu[i] ~ dnorm(0,0.0001) # random effect on baselines
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

	# Model for baseline effects ###################################
	# Adapted from Program 1 of NICE DSU TSD 5
	for (i in 1:ns.base){ # LOOP THROUGH STUDIES
		r.base[i] ~ dbin(p.base[i],n.base[i]) # Likelihood
		logit(p.base[i]) <- mu.base[i] # Log-odds of response
		mu.base[i] ~ dnorm(m,tau.m) # Random effects model
	}
	m ~ dnorm(0,.0001) # vague prior for mean
	var.m <- 1/tau.m # between-trial variance
	tau.m <- pow(sd.m,-2) # between-trial precision = (1/between-trial variance)
	sd.m ~ dunif(0,5) # vague prior for between-trial SD


	# Model for single-arm studies ###################################
	# Prevent feedback to baseline model (these two may not be necessary as mu.single is cut below
	m.cut<-cut(m)
	tau.m.cut<-cut(tau.m)
	for(i in 1:ns.single)
	{
		mu.single[i]~dnorm(m.cut,tau.m.cut) # Sampled from baseline random effects model	
		mu.single.cut[i]<-cut(mu.single[i]) # Prevent feedback to baseline model
		r.single[i]~dbin(p.single[i],n.single[i])
		logit(p.single[i])<-mu.single.cut[i]+delta.single[i]
		delta.single[i]<-d[t.single[i]] # Treatment effect relative to reference

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
	for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects 
}



