# Simulation study to investigate ALM and reference prediction
# Script to define parameters of the 6 scenarios
# Howard Thom 21-December-2019

# v2 doubles number of scenarios. One set has gamma being Normal(mean=0.5, sd=1) the other has gamma equal to 0

scenario.params<-list()

n.scenarios <- 12

# Scenarios 1 and 2
scenario.params[[1]]<-list()
scenario.params[[1]]$nt<-5
scenario.params[[1]]$ns<-5
scenario.params[[1]]$na<-c(2,2,2,2,3)
scenario.params[[1]]$tr<-t(matrix(c(1,2,NA,1,2,NA,
						1,3,NA,1,3,NA,
						1,2,3),nrow=3))
scenario.params[[1]]$ns.base<-5
scenario.params[[1]]$t.single<-c(4,4,4,4,4,5,5,5,5,5)
scenario.params[[1]]$ns.single<-10
scenario.params[[1]]$ns.disc<-5
scenario.params[[1]]$na.disc<-c(2,2,2,2,2)
scenario.params[[1]]$t.disc<-t(matrix(c(4,5,NA,4,5,NA,
						4,5,NA,4,5,NA,
						4,5,NA),nrow=3))

# Only the connected portion changes
scenario.params[[2]] <- scenario.params[[3]] <- scenario.params[[4]] <- 
	scenario.params[[5]] <- scenario.params[[6]] <- scenario.params[[1]]

# Scenarios 3 and 4
scenario.params[[3]]$ns <- 15
scenario.params[[3]]$na <- c(rep(2, 5),rep(2, 5), rep(3, 5))
scenario.params[[3]]$tr <- t(matrix(c(rep(c(1, 2, NA), 5),
						rep(c(1, 3, NA), 5),
						rep(c(1, 2, 3), 5)), nrow=3))
scenario.params[[4]] <- scenario.params[[3]]

# Scenarios 5 and 6
scenario.params[[5]]$ns <- 50
scenario.params[[5]]$na <- c(rep(2, 20),rep(2, 20), rep(3, 10))
scenario.params[[5]]$tr <- t(matrix(c(rep(c(1, 2, NA), 20),
						rep(c(1, 3, NA), 20),
						rep(c(1, 2, 3), 10)), nrow=3))
scenario.params[[6]] <- scenario.params[[5]]


# Strength of covariate
# Odd scenarios are weak
for(i.scenario in 2*c(0:2)+1){
	scenario.params[[i.scenario]]$beta.mean <- 0.1
	scenario.params[[i.scenario]]$beta.sd <- 1
}
# Even scenarios are strong
for(i.scenario in 2*c(1:3)){
	scenario.params[[i.scenario]]$beta.mean <- 1
	scenario.params[[i.scenario]]$beta.sd <- 1
}

# Now add the gamma scenarios
# Scenarios with a random difference in baseline response
for(i.scenario in 1:6) {
  scenario.params[[i.scenario]]$gamma.mean <- 0.5
  scenario.params[[i.scenario]]$gamma.sd <- 1
}
# Scenarios with no difference in baselien response
for(i.scenario in 7:12) {
  scenario.params[[i.scenario]] <- scenario.params[[i.scenario-6]]
  scenario.params[[i.scenario]]$gamma.mean <- 0
  scenario.params[[i.scenario]]$gamma.sd <- 0
}


