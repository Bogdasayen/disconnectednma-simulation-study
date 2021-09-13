# Simulation study to investigate ALM and reference prediction
# Function to generate artificial datasets for random effects on baseline
# and ALM simulation study
# Howard Thom 5-November-2019

# V3 allows gamma to vary

# Data used by single-arm studies analyses
# Data are same as TSD format: ns, nt, na, r, n, t, x
# The data on baseline arms (from RCTs) are
# ns.base : number of baseline arms
# r.base : number of events in baseline arms
# n.base : number of patients in baseline arms
# x.base : covariate value for baseline arms
# x.base.mean : average covariate value for baseline arms
# The data on single arm studies are:
# ns.single : number of single arm studies
# r.single : number of events in single arm studies
# n.single : number of patients in single arm studies
# x.single : covariate value for single arm studies
# t.single : treatment in single arm study


# These are likely not going to be used.
# We'll use the same simulated data for fixed and random effects analyses
# Extra data for random effects models is added during model estimation
# sd.connected.mean and sd.connected.tau are the mean and precision of the sd in the connected components.



# Logistic link function
logit<-function(x)
{
	return(log(x/(1-x)))
}

# Inverse of logit
expit<-function(x)
{
	return(1/(1+exp(-x)))
}

# Generate the connected portion of the network
simulate.bugs.data<-function(i.scenario=1)
{
	# Global simulation parameters
	# These define the structure of RCTs and single-arm studies
	# that are shared across simulations
	# For all treatments
	nt <- scenario.params[[i.scenario]]$nt
	
	# For the connected RCTs
	ns <- scenario.params[[i.scenario]]$ns
	na <- scenario.params[[i.scenario]]$na
	tr <- scenario.params[[i.scenario]]$tr
	beta.mean <- scenario.params[[i.scenario]]$beta.mean
	beta.sd <- scenario.params[[i.scenario]]$beta.sd

	# For single arm studies
	ns.base<-scenario.params[[i.scenario]]$ns.base
	ns.single<-scenario.params[[i.scenario]]$ns.single
	t.single<-scenario.params[[i.scenario]]$t.single
	gamma.mean <- scenario.params[[i.scenario]]$gamma.mean
	gamma.sd <- scenario.params[[i.scenario]]$gamma.sd

	# For disconnected RCTs
	ns.base<-scenario.params[[i.scenario]]$ns.base
	na.disc<-scenario.params[[i.scenario]]$na.disc
	t.disc<-scenario.params[[i.scenario]]$t.disc
	ns.disc<-scenario.params[[i.scenario]]$ns.disc


	# Model parameters
	m <- rnorm(1, mean = 0.5, sd = 1)
	gamma <- rnorm(1, mean = gamma.mean, sd = gamma.sd)
	# Covariate effect (scenario dependent)
	beta <- rnorm(1, mean = beta.mean, sd = beta.sd)
	# Covariate for each arm of connected RCTs
	x <- matrix(rnorm(ns * max(na), mean = 0.5, sd = 1), nrow = ns)
	x[na != 3, 3] <- NA
	# Covariate for each arm of disconnected RCTs
	x.disc <- matrix(rnorm(ns.disc * max(na.disc), mean = 0.5, sd = 1), nrow = ns.disc)
	# Covariate for each of the single arm studies
	x.single <- rnorm(ns.single, mean = 0.5, sd =1)
	# Treatment effects
	d <- rnorm(nt, mean = 0, sd = 1)
	# Baseline corresponds to treatment 1
	d[1] <- 0

	# Study level baseline effect for connected RCTs
	mu <- m + x[, 1] * beta + rnorm(ns, mean = 0, sd = 1) 
	# Study level numbers of patients fixed at n.patients 
	# Note that this code only works if maximum number of arms is 3
	n <- matrix(n.patients, nrow = ns, ncol = max(na))
	n[na != 3, 3] <- NA
	# Study level probabilities and numbers of event
	p <- r <- matrix(NA, nrow = ns, ncol = max(na))
	for(i in 1:ns)
	{
		for(j in 1:na[i])
		{
		p[i, j] <- expit(mu[i] + d[tr[i, j]])
		r[i, j] <- rbinom(1, size = n[i, j], prob = p[i, j])
		}
	}
	# Do a continuity correction on r
	r[r==0 & !is.na(r)]<-1

	# Extract the baseline arm data from the connected RCTs
	ns.base <- sum(tr == 1, na.rm = TRUE)
	r.base <- r[tr==1 & !is.na(tr)]
	n.base <- n[tr==1 & !is.na(tr)]
	x.base <- x[tr==1 & !is.na(tr)]
	x.base.mean <- mean(x.base)

	# Study level baseline effect for disconnected RCTs
	mu.disc <- m + x.disc[, 1] * beta + x.disc[, 1] * gamma + rnorm(ns.disc, mean = 0, sd = 1)
	# Note this code assumes all disconnected studies have only 2 arms
	n.disc <- matrix(100, nrow = ns.disc, ncol = max(na.disc))
	# Study level probabilities and numbers of event
	p.disc <- r.disc <- matrix(NA, nrow = ns.disc, ncol = max(na.disc))
	for(i in 1:ns.disc)
	{
		for(j in 1:na.disc[i])
		{
		p.disc[i, j] <- expit(mu.disc[i] + d[t.disc[i, j]])
		r.disc[i, j] <- rbinom(1, size = n.disc[i, j], prob = p.disc[i, j])
		}
	}
	# Do a continuity correction on r
	r.disc[r.disc==0 & !is.na(r.disc)]<-1

	# Study level baseline effect for disconnected RCTs
	mu.single <- m + x.single * beta + x.single * gamma + rnorm(ns.single, mean = 0, sd = 1)
	n.single <- rep(100, ns.single)
	p.single <- r.single <- rep(NA, ns.single)
	for(i in 1:ns.single)
	{
		r.single[i] <- rbinom(1, size = n.single[i], prob = expit(mu.single[i] + d[t.single[i]]))
	}
	# Do a continuity correction on r
	r.single[r.single==0 & !is.na(r.single)]<-1


	# Create the dataset for OpenBUGS
	b.data <- list(ns = ns, nt = nt, r = r, n = n, t = tr, x = x, na = na,
			ns.base = ns.base, r.base = r.base, n.base = n.base, x.base = x.base,
			x.base.mean = x.base.mean,
			ns.single = ns.single, r.single = r.single, n.single = n.single, t.single = t.single, x.single = x.single,
			ns.disc = ns.disc, r.disc = r.disc, n.disc = n.disc, t.disc = t.disc, x.disc = x.disc, na.disc = na.disc)

	# Store the simulation parameters (only the d will be compared to estimated values though)
	sim.params <- list(d = d, m = m, beta = beta, gamma = gamma)

	# Export the BUGS data and simulation parameters
	return(list(sim.params = sim.params,
			b.data = b.data))
}




