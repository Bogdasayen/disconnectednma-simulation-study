# Baseline prediction models with disconnected RCTs
# Howard Thom 6-October-2017
# Found that N(0,0.298) and Unif(0,2) priors helped achieve convergence.

# Changes from previous versions
# Uses only reference treatment arms to fit the baseline model
# Uses 'cut()' to prevent feedback from single-arm studies to baseline model
# Naively pools disconnected RCT evidence on treatment effects (not sure enough evidence to fit hierarchical model)

# Change from version 1
# Ensure each arm k of a disconnected RCT i informs d[t[i,k]]-t[1]
# instead of d[t[i,k]]-d[t[i,1]].
# The reason is that the mu[i] represents treatment t[1] (Reference), not t[i,1] (baseline)
# Note slight difference in multi-arm correction as now delta[t[i,1]] is not zero

# V3 puts the random effect around a fixed intercept with regression
# v3 also centres the covariates at x.base.mean for both fixed and random effects models


# Data are same as TSD format: ns, nt, na, r, n, t
# The data on baseline arms (from RCTs) are
# ns.base : number of baseline arms
# r.base : number of events in baseline arms
# n.base : number of patients in baseline arms
# x.base : Matrix of covariates with one row per covariate
# The data on disconnected networks are (as in standard TSD):
# ns.disc, nt.disc, na.disc, r.disc, n.disc, t.disc
# x.disc is matrix of covariates for disconnected RCTs.
# cov.index is the index of the covariates to use in x.disc and x.base

# NOTE: All x.base and x.disc must be defined. Suggest using mean covariate value when not reported.
# NOTE: Naively pools disconnected and connected RCT evidence on treatment effects (not sure enough evidence to fit hierarchical model)

# Model with random effecs on the baseline
# Binomial likelihood, logit link 
# Simultaneous baseline and treat effects model for multi-arm trials 
# Includes disconnected RCTs
model.random.effects.baseline.disc.1cov.re<-function()
{ # *** PROGRAM STARTS 
	# Model for RCTs ##############################################
	for(i in 1:ns){ # LOOP THROUGH STUDIES 
		w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm 
		delta[i,1] <- 0 # treatment effect is zero for control arm 
		mu[i] ~ dnorm(0,0.298) # Baseline are a nuisance for RCT evidence
		for (k in 1:na[i]) { # LOOP THROUGH ARMS 
			r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood 
			logit(p[i,k]) <- mu[i] + delta[i,k] # model for linear predictor 
			rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
			dev.NA[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution including NAs 
			+ (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k]))) 
			dev[i,k] <- dev.NA[i,k]*(1-equals(n[i,1],1)) # Deviance contribution with correction for NAs 
		} 
		resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial 
		for (k in 2:na[i]) { # LOOP THROUGH ARMS 
			delta[i,k] ~ dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions 
			md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction) 
			taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction) 
			w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs 
			sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials 
		} 
	}

	# Model for baseline effects ###################################
	# Adapted from Program 1 of NICE DSU TSD 5
	for (i in 1:ns.base){ # LOOP THROUGH STUDIES
		r.base[i] ~ dbin(p.base[i],n.base[i]) # Likelihood
		logit(p.base[i]) <- mu.base[i] 
		mu.base[i] ~ dnorm(mu.base.mean[i],tau.m) # Random effects model
		mu.base.mean[i] <- m + (x.base[i]-x.base.mean)*beta.base[1] # Prediction of mean effect # Log-odds of response
	}
	beta.base~dnorm(0,0.298) #  vague prior for covariate effects

	m ~ dnorm(0,.298) # vague prior for mean
	var.m <- 1/tau.m # between-trial variance
	tau.m <- pow(sd.m,-2) # between-trial precision = (1/between-trial variance)
	sd.m ~ dunif(0,2) # vague prior for between-trial SD


	# Model for disconnected RCTs #######################################
	# Prevent feedback to baseline model (these two may not be necessary as mu.disc is cut below
	m.cut<-cut(m)
	tau.m.cut<-cut(tau.m)
	beta.base.cut<-cut(beta.base)
	for(i in 1:ns.disc){ # LOOP THROUGH STUDIES 
		mu.disc[i]~dnorm(m.cut,tau.m.cut) # Sampled from baseline random effects model
		mu.disc.cut[i]<-cut(mu.disc[i]) # Prevent feedback to baseline model
		for (k in 1:na.disc[i]) { # LOOP THROUGH ARMS 
			r.disc[i,k] ~ dbin(p.disc[i,k],n.disc[i,k]) # binomial likelihood 
			logit(p.disc[i,k]) <- mu.disc[i, k]
			mu.disc[i, k] ~ dnorm(mu.disc.mean[i, k], tau.m.cut)
			mu.disc.mean[i, k] <- m.cut + delta.disc[i,k] + (x.disc[i,k]-x.base.mean)*beta.base.cut # model for linear predictor 
			rhat.disc[i,k] <- p.disc[i,k] * n.disc[i,k] # expected value of the numerators 
			dev.NA.disc[i,k] <- 2 * (r.disc[i,k] * (log(r.disc[i,k])-log(rhat.disc[i,k])) #Deviance contribution including NAs 
			+ (n.disc[i,k]-r.disc[i,k]) * (log(n.disc[i,k]-r.disc[i,k]) - log(n.disc[i,k]-rhat.disc[i,k]))) 
			dev.disc[i,k] <- dev.NA.disc[i,k]*(1-equals(n.disc[i,1],1)) # Deviance contribution with correction for NAs 
		} 
		resdev.disc[i] <- sum(dev.disc[i,1:na.disc[i]]) # summed residual deviance contribution for this trial 

		delta.disc[i,1] ~ dnorm(d[t.disc[i,1]], tau.disc) # Handle k=1 separately as sum dissappears.
		w.disc[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm 

		for (k in 2:na.disc[i]) { # LOOP THROUGH ARMS (including control arm)
			delta.disc[i,k] ~ dnorm(md.disc[i,k],taud.disc[i,k]) # trial-specific LOR distributions 
			md.disc[i,k] <- d[t.disc[i,k]] + sw.disc[i,k] # mean of LOR distributions (with multi-arm trial correction) 
			taud.disc[i,k] <- tau.disc *2*k/(k+1) # precision of LOR distributions (with multi-arm trial correction) 
			w.disc[i,k] <- (delta.disc[i,k] - d[t.disc[i,k]] ) # adjustment for multi-arm RCTs 
			sw.disc[i,k] <- sum(w.disc[i,1:(k-1)])/k # cumulative adjustment for multi-arm trials 
		} 
	}

	# Calculate deviance ##############################################
	totresdev.disc<-sum(resdev.disc[]) # Total residual deviance for disconnected RCTs
	totresdev.rct <- sum(resdev[]) # Total residual deviance for RCTs
	totresdev<-totresdev.disc+totresdev.rct #Total Residual Deviance 
	
	# Specify remaining priors
	d[1]<-0 # treatment effect is zero for reference treatment 
	for (k in 2:nt){ d[k] ~ dnorm(0,.298) } # vague priors for treatment effects 
	sd ~ dunif(0,2) # vague prior for between-trial SD 
	tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance) 

	sd.disc ~ dunif(0,2) # vague prior for between-trial SD 
	tau.disc <- pow(sd.disc,-2) # between-trial precision = (1/between-trial variance) 
}


# The random effects on baseline and fixed treatment effect model, including disconnected RCTs
model.random.effects.baseline.disc.1cov.fe<-function()
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

	# Model for baseline effects ###################################
	# Adapted from Program 1 of NICE DSU TSD 5
	for (i in 1:ns.base){ # LOOP THROUGH STUDIES
		r.base[i] ~ dbin(p.base[i],n.base[i]) # Likelihood
		logit(p.base[i]) <- mu.base[i] 
		mu.base[i] ~ dnorm(mu.base.mean[i],tau.m) # Random effects model
		mu.base.mean[i] <- m + (x.base[i]-x.base.mean)*beta.base # Prediction of mean effect # Log-odds of response
	}
	beta.base~dnorm(0,0.01) #0.298 #  vague prior for covariate effects
	m ~ dnorm(0,.01) # vague prior for mean
	var.m <- 1/tau.m # between-trial variance
	tau.m <- pow(sd.m,-2) # between-trial precision = (1/between-trial variance)
	sd.m ~ dunif(0,2) # vague prior for between-trial SD


	# Model for disconnected RCTs ###################################
	# Prevent feedback to baseline model (these two may not be necessary as mu.single is cut below
	m.cut<-cut(m)
	tau.m.cut<-cut(tau.m)
	beta.base.cut<-cut(beta.base)

	for(i in 1:ns.disc){ # LOOP THROUGH STUDIES 
		for (k in 1:na.disc[i]) { # LOOP THROUGH ARMS 
			r.disc[i,k] ~ dbin(p.disc[i,k],n.disc[i,k]) # binomial likelihood 
			logit(p.disc[i,k]) <- mu.disc[i,k]
			mu.disc[i,k] ~ dnorm(mu.disc.mean[i,k], tau.m.cut)
			mu.disc.mean[i,k] <- m.cut + delta.disc[i,k] + (x.disc[i,k]-x.base.mean)*beta.base.cut #+ x.disc[cov.index[2],i,k]*beta.base.cut[2] #+ x.disc[cov.index[3],i,k]*beta.base.cut[3]  # model for linear predictor 
			delta.disc[i,k]<-d[t.disc[i,k]] - d[t.disc[i,1]] # model for linear predictor 
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
	for (k in 2:nt){ d[k] ~ dnorm(0,.01) } # vague priors for treatment effects 
}



