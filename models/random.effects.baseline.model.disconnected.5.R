# Random effects on baseline models with disconnected RCTs
# Howard Thom 6-October-2017

# Changes from previous versions
# Uses only reference treatment arms to fit the baseline model
# Uses 'cut()' to prevent feedback from single-arm studies to baseline model
# Naively pools disconnected RCT evidence on treatment effects (not sure enough evidence to fit hierarchical model)

# Change from version 1
# Ensure each arm k of a disconnected RCT i informs d[t[i,k]]-d[1]
# instead of d[t[i,k]]-d[t[i,1]].
# The reason is that the mu[i] represents treatment t[1] (Reference), not t[i,1] (baseline)
# Note slight difference in multi-arm correction as now delta[t[i,1]] is not zero

# V2 changes priors on mu to dnorm(0,0.01) rather than dnorm(0,0.298). Vaguer priors do not converge.
# V4 uses separate standard deviation for random effecs models in connected and single-arm portions
# V5 removes the random effects code as this has to be put in a separate txt file where the dnorm(,)I(,) syntax can be used.

# Data are same as TSD format: ns, nt, na, r, n, t
# The data on baseline arms (from RCTs) are
# ns.base : number of baseline arms
# r.base : number of events in baseline arms
# n.base : number of patients in baseline arms
# The data on disconnected networks are (as in standard TSD):
# ns.disc, nt.disc, na.disc, r.disc, n.disc, t.disc

# sd.connected.mean and sd.connected.tau are the mean and precision of the sd in the connected components.
# These are used as informative priors on sd.disc in the random effects models

# The random effects on baseline and fixed treatment effect model, including disconnected RCTs
model.random.effects.baseline.disc.fe<-function()
{
	# Model for RCTs ###################################
	for(i in 1:ns){ # LOOP THROUGH STUDIES 
		mu[i] ~ dnorm(0,0.001) # vague prior on independent baselines
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
	m ~ dnorm(0,.001) # vague prior for mean
	var.m <- 1/tau.m # between-trial variance
	tau.m <- pow(sd.m,-2) # between-trial precision = (1/between-trial variance)
	sd.m ~ dunif(0,5) # vague prior for between-trial SD


	# Model for disconnected RCTs ###################################
	# Prevent feedback to baseline model (these two may not be necessary as mu.single is cut below
	m.cut<-cut(m)
	tau.m.cut<-cut(tau.m)
	for(i in 1:ns.disc){ # LOOP THROUGH STUDIES 
		mu.disc[i]~dnorm(m.cut,tau.m.cut) # Sampled from baseline random effects model	
		mu.disc.cut[i]<-cut(mu.disc[i]) # Prevent feedback to baseline model
		for (k in 1:na.disc[i]) { # LOOP THROUGH ARMS 
			r.disc[i,k] ~ dbin(p.disc[i,k],n.disc[i,k]) # Binomial likelihood 
			logit(p.disc[i,k]) <- mu.disc.cut[i] + delta.disc[i,k]
			delta.disc[i,k]<-d[t.disc[i,k]] - d[1] # model for linear predictor 
			rhat.disc[i,k] <- p.disc[i,k] * n.disc[i,k] # expected value of the numerators 
			dev.disc[i,k] <- 2 * (r.disc[i,k] * (log(r.disc[i,k])-log(rhat.disc[i,k])) 
			+ (n.disc[i,k]-r.disc[i,k]) * (log(n.disc[i,k]-r.disc[i,k]) - log(n.disc[i,k]-rhat.disc[i,k]))) #Deviance contribution 
		} 
	resdev.disc[i] <- sum(dev.disc[i,1:na.disc[i]]) # summed residual deviance contribution for this trial 
	} 

	# Calculate deviance ###################################
	totresdev.disc<-sum(resdev.disc[]) # Total residual deviance for disconnected RCTs
	totresdev.rct <- sum(resdev[]) # Total residual deviance for RCTs
	totresdev<-totresdev.disc+totresdev.rct #Total Residual Deviance 
	
	# Priors for remaining parameters
	d[1]<-0 # treatment effect is zero for reference treatment 
	for (k in 2:nt){ d[k] ~ dnorm(0,.001) } # vague priors for treatment effects 
}



