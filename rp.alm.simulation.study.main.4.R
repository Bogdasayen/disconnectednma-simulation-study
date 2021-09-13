# Simulation study to investigate ALM and reference prediction
# Main script to generate datasets, run analyses, and estimate bias and coverage
# Howard Thom 21-December-2019


# V2 added the gamma scenarios
# V3 includes random effects
# v4 includes plots and allows patient numbers to vary

# TODO Extract results from random effects simulations


library(openxlsx)

source("define.simulation.parameters.2.R")
source("analyse.simulation.data.4.R") 
source("generate.simulation.data.3.R")
source("plot.simulation.results.1.R")

n.simulations <- 100

# OpenBUGS settings
n.chains <- 2 # The simulation code is not set up for anything other than 2
num.sims <- 1000#0 * n.chains
burn.in <- 1000#0 * n.chains  (recommend 30000 for random effects, plus 10000 samples and n.thin=2)

n.patients <- 100

# Choose whether to do fixed and/or random effects
random.effects <- TRUE
fixed.effects <- TRUE

bias.tables <- bias.tables.re <- list()
coverage.tables <- coverage.tables.re <- list()
run.time<-system.time({
for(i.scenario in 1:n.scenarios)
{
  
	bias.tables[[i.scenario]] <- bias.tables.re[[i.scenario]] <-
		coverage.tables[[i.scenario]] <- coverage.tables.re[[i.scenario]] <-
		matrix(NA, nrow = n.simulations, ncol = 9, 
		dimnames = list(NULL, c("RCT Only connected", "ALM single connected", "ALM disc connected",
			"RP single connected", "RP disc connected", 
			"ALM single disconnected", "ALM disc disconnected",
			"RP single disconnected", "RP disc disconnected")))

	for(i.simulation in 1:n.simulations)
	{
		print(paste("Scenario",i.scenario,", simulation",i.simulation))
		# Simulate the disconnected networks for NMA
		temp.data <- simulate.bugs.data(i.scenario)
		
		# Use OpenBUGS to analyse the simulated data
		if(fixed.effects) {
		  print("Fixed effects")
		  temp.results <- analyse.simulation.data(b.data = temp.data$b.data,
		                                          num.sims = num.sims, burn.in = burn.in, do.debug = FALSE, random.effects = FALSE)  
		  
		  # Bias is the difference between mean estimates and truth
		  bias.tables[[i.scenario]][i.simulation, c("RCT Only connected", "ALM single connected", "ALM disc connected",
		                                            "RP single connected", "RP disc connected")] <- 
		    rowMeans(temp.results$mean.table[, 1:2] - t(matrix(rep(temp.data$sim.params$d[2:3], 5), nrow = 2)))
		  
		  bias.tables[[i.scenario]][i.simulation, c("ALM single disconnected", "ALM disc disconnected",
		                                            "RP single disconnected", "RP disc disconnected")] <- 
		    rowMeans(temp.results$mean.table[-1 , 3:4] - t(matrix(rep(temp.data$sim.params$d[4:5], 4), nrow = 2)))
		  
		  # Coverage is logical whether or not the 95% CrI includes truth
		  # Bias is the difference between mean estimates and truth
		  coverage.tables[[i.scenario]][i.simulation, c("RCT Only connected", "ALM single connected", "ALM disc connected",
		                                                "RP single connected", "RP disc connected")] <- 
		    rowMeans((temp.results$ll.table[, 1:2] - t(matrix(rep(temp.data$sim.params$d[2:3], 5), nrow = 2)))<=0 &
		               (temp.results$ul.table[, 1:2] - t(matrix(rep(temp.data$sim.params$d[2:3], 5), nrow = 2)))>=0)
		  
		  coverage.tables[[i.scenario]][i.simulation, c("ALM single disconnected", "ALM disc disconnected",
		                                                "RP single disconnected", "RP disc disconnected")] <- 
		    rowMeans((temp.results$ll.table[-1, 3:4] - t(matrix(rep(temp.data$sim.params$d[4:5], 4), nrow = 2)))<=0 &
		               (temp.results$ul.table[-1, 3:4] - t(matrix(rep(temp.data$sim.params$d[4:5], 4), nrow = 2)))>=0)
		  
		} # End fixed effects
		
		if(random.effects) {
		  print("Random effects")
		  temp.results.re <- analyse.simulation.data(b.data = temp.data$b.data,
		                                             num.sims = num.sims, burn.in = burn.in, do.debug = FALSE, random.effects = TRUE)
		  
		  # Bias is the difference between mean estimates and truth
		  bias.tables.re[[i.scenario]][i.simulation, c("RCT Only connected", "ALM single connected", "ALM disc connected",
		                                            "RP single connected", "RP disc connected")] <- 
		    rowMeans(temp.results.re$mean.table[, 1:2] - t(matrix(rep(temp.data$sim.params$d[2:3], 5), nrow = 2)))
		  
		  bias.tables.re[[i.scenario]][i.simulation, c("ALM single disconnected", "ALM disc disconnected",
		                                            "RP single disconnected", "RP disc disconnected")] <- 
		    rowMeans(temp.results.re$mean.table[-1 , 3:4] - t(matrix(rep(temp.data$sim.params$d[4:5], 4), nrow = 2)))
		  
		  # Coverage is logical whether or not the 95% CrI includes truth
		  # Bias is the difference between mean estimates and truth
		  coverage.tables.re[[i.scenario]][i.simulation, c("RCT Only connected", "ALM single connected", "ALM disc connected",
		                                                "RP single connected", "RP disc connected")] <- 
		    rowMeans((temp.results.re$ll.table[, 1:2] - t(matrix(rep(temp.data$sim.params$d[2:3], 5), nrow = 2)))<=0 &
		               (temp.results.re$ul.table[, 1:2] - t(matrix(rep(temp.data$sim.params$d[2:3], 5), nrow = 2)))>=0)
		  
		  coverage.tables.re[[i.scenario]][i.simulation, c("ALM single disconnected", "ALM disc disconnected",
		                                                "RP single disconnected", "RP disc disconnected")] <- 
		    rowMeans((temp.results.re$ll.table[-1, 3:4] - t(matrix(rep(temp.data$sim.params$d[4:5], 4), nrow = 2)))<=0 &
		               (temp.results.re$ul.table[-1, 3:4] - t(matrix(rep(temp.data$sim.params$d[4:5], 4), nrow = 2)))>=0)
		  
		  
		} # End random effects
	
	}
  
}
})# End system time

# If not running the simulation afresh, can use results of previous run
load(file=paste0(baseline.directory,"/code/simulation.results.",n.simulations,".v2.rda"))
load(file=paste0(baseline.directory,"/code/simulation.results.re.",n.simulations,".v2.rda"))

# Summarise each scenario
if(fixed.effects) {
  # Construct matrices to store all scenarios and both single and disconnected network results
  connected.summary.matrix.names <- list(c("RCT only","ALM single", "ALM disc", "RP single", "RP disc"), paste(rep("scenario",n.scenarios),c(1:n.scenarios)))
  disconnected.summary.matrix.names <- list(c("ALM single", "ALM disc", "RP single", "RP disc"), paste(rep("scenario",n.scenarios),c(1:n.scenarios)))
  bias.summary.connected <- bias.mean.connected <- bias.ll.connected <- bias.ul.connected <- matrix(nrow = 5, ncol = n.scenarios, dimnames = connected.summary.matrix.names )
  coverage.summary.connected <- coverage.mean.connected <- coverage.ll.connected <- coverage.ul.connected <- matrix(nrow = 5, ncol = n.scenarios, dimnames = connected.summary.matrix.names )
  bias.summary.disconnected <- bias.mean.disconnected <- bias.ll.disconnected <- bias.ul.disconnected <- matrix(nrow = 4, ncol = n.scenarios, dimnames = disconnected.summary.matrix.names )
  coverage.summary.disconnected <- coverage.mean.disconnected <- coverage.ll.disconnected <- coverage.ul.disconnected  <- matrix(nrow = 4, ncol = n.scenarios, dimnames = disconnected.summary.matrix.names )
  
  for(i.scenario in 1:n.scenarios)
  {
    # Mean and 95% limits for use in figure
    bias.mean.connected[, i.scenario] <- colMeans(bias.tables[[i.scenario]], na.rm = TRUE)[1:5]
    bias.mean.disconnected[, i.scenario] <- colMeans(bias.tables[[i.scenario]], na.rm = TRUE)[6:9]
    coverage.mean.connected[, i.scenario] <- colMeans(coverage.tables[[i.scenario]], na.rm = TRUE)[1:5]
    coverage.mean.disconnected[, i.scenario] <- colMeans(coverage.tables[[i.scenario]], na.rm = TRUE)[6:9]
    
    bias.ll.connected[, i.scenario] <- apply(bias.tables[[i.scenario]], c(2), mean)[1:5] - 1.96 * apply(bias.tables[[i.scenario]], c(2), sd)[1:5]/sqrt(n.simulations)
    bias.ll.disconnected[, i.scenario] <- apply(bias.tables[[i.scenario]], c(2), mean)[6:9] - 1.96 * apply(bias.tables[[i.scenario]], c(2), sd)[6:9]/sqrt(n.simulations)
    coverage.ll.connected[, i.scenario] <-  apply(coverage.tables[[i.scenario]], c(2), mean)[1:5] - 1.96 * apply(coverage.tables[[i.scenario]], c(2), sd)[1:5]/sqrt(n.simulations)
    coverage.ll.disconnected[, i.scenario] <-  apply(coverage.tables[[i.scenario]], c(2), mean)[6:9] - 1.96 * apply(coverage.tables[[i.scenario]], c(2), sd)[6:9]/sqrt(n.simulations)
    
    bias.ul.connected[, i.scenario] <- apply(bias.tables[[i.scenario]], c(2), mean)[1:5] + 1.96 * apply(bias.tables[[i.scenario]], c(2), sd)[1:5]/sqrt(n.simulations)
    bias.ul.disconnected[, i.scenario] <- apply(bias.tables[[i.scenario]], c(2), mean)[6:9] + 1.96 * apply(bias.tables[[i.scenario]], c(2), sd)[6:9]/sqrt(n.simulations)
    coverage.ul.connected[, i.scenario] <-  apply(coverage.tables[[i.scenario]], c(2), mean)[1:5] + 1.96 * apply(coverage.tables[[i.scenario]], c(2), sd)[1:5]/sqrt(n.simulations)
    coverage.ul.disconnected[, i.scenario] <-  apply(coverage.tables[[i.scenario]], c(2), mean)[6:9] + 1.96 * apply(coverage.tables[[i.scenario]], c(2), sd)[6:9]/sqrt(n.simulations)
    
    
  }
  # Ensure coverage does not extend beyond 0, 1
  coverage.ll.connected[coverage.ll.connected < 0] <- 0
  coverage.ll.disconnected[coverage.ll.disconnected < 0] <- 0
  coverage.ul.connected[coverage.ul.connected > 1] <- 1
  coverage.ul.disconnected[coverage.ul.disconnected > 1] <- 1
  for(i.scenario in 1:n.scenarios) 
  {
    # Summary for use in tables
    bias.summary.connected[, i.scenario] <- paste0(format(bias.mean.connected[, i.scenario], digits = 2), " (",
                                                   format(bias.ll.connected[, i.scenario], digits = 2), ", ",
                                                   format(bias.ul.connected[, i.scenario], digits = 2), ")")
    
    bias.summary.disconnected[, i.scenario] <- paste0(format(bias.mean.disconnected[, i.scenario], digits = 2), " (",
                                                      format(bias.ll.disconnected[, i.scenario], digits = 2), ", ",
                                                      format(bias.ul.disconnected[, i.scenario], digits = 2), ")")
    coverage.summary.connected[, i.scenario] <- paste0(format(coverage.mean.connected[, i.scenario], digits = 2), " (",
                                                       format(coverage.ll.connected[, i.scenario], digits = 2), ", ",
                                                       format(coverage.ul.connected[, i.scenario], digits = 2), ")")
    coverage.summary.disconnected[, i.scenario] <- paste0(format(coverage.mean.disconnected[, i.scenario], digits = 2), " (",
                                                          format(coverage.ll.disconnected[, i.scenario], digits = 2), ", ",
                                                          format(coverage.ul.disconnected[, i.scenario], digits = 2), ")")
    
  }
  
  plot.simulation.results(bias.mean.connected,bias.mean.disconnected,bias.ll.connected, bias.ll.disconnected, bias.ul.connected, bias.ul.disconnected,
                               coverage.mean.connected,coverage.mean.disconnected,coverage.ll.connected, coverage.ll.disconnected, coverage.ul.connected, coverage.ul.disconnected,
                          start.of.filename = paste0(baseline.directory,"/results/plots/simulation.fe.",n.simulations, ".", n.patients, "patients"))
  
  save(bias.summary.connected,bias.summary.disconnected,coverage.summary.connected,coverage.summary.disconnected,bias.tables,coverage.tables,
       bias.mean.connected,bias.mean.disconnected,bias.ll.connected, bias.ll.disconnected, bias.ul.connected, bias.ul.disconnected,
       coverage.mean.connected,coverage.mean.disconnected,coverage.ll.connected, coverage.ll.disconnected, coverage.ul.connected, coverage.ul.disconnected,
       file=paste0(baseline.directory,"/code/simulation.results.fe.",n.simulations,".", n.patients, "patients.v3.rda"))
  
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "Bias conn")
  writeData(wb, "Bias conn", cbind(rownames(bias.summary.connected), bias.summary.connected), startRow = 1, startCol = 1)
  addWorksheet(wb, "Bias disc")
  writeData(wb, "Bias disc", cbind(rownames(bias.summary.disconnected), bias.summary.disconnected), startRow = 1, startCol = 1)
  addWorksheet(wb, "Coverage conn")
  writeData(wb, "Coverage conn", cbind(rownames(coverage.summary.connected), coverage.summary.connected), startRow = 1, startCol = 1)
  addWorksheet(wb, "Coverage disc")
  writeData(wb, "Coverage disc", coverage.summary.disconnected, startRow = 1, startCol = 2)
  saveWorkbook(wb, file = paste0(baseline.directory,"/results/simulation.results.fe.",n.simulations,".", n.patients, "patients.v3.xlsx"), overwrite = TRUE)    

}

if(random.effects) {
  # Construct matrices to store all scenarios and both single and disconnected network results
  connected.summary.matrix.names <- list(c("RCT only","ALM single", "ALM disc", "RP single", "RP disc"), paste(rep("scenario",n.scenarios),c(1:n.scenarios)))
  disconnected.summary.matrix.names <- list(c("ALM single", "ALM disc", "RP single", "RP disc"), paste(rep("scenario",n.scenarios),c(1:n.scenarios)))
  bias.summary.connected.re <- bias.mean.connected.re <- bias.ll.connected.re <- bias.ul.connected.re <- matrix(nrow = 5, ncol = n.scenarios, dimnames = connected.summary.matrix.names )
  coverage.summary.connected.re <- coverage.mean.connected.re <- coverage.ll.connected.re <- coverage.ul.connected.re <- matrix(nrow = 5, ncol = n.scenarios, dimnames = connected.summary.matrix.names )
  bias.summary.disconnected.re <- bias.mean.disconnected.re <- bias.ll.disconnected.re <- bias.ul.disconnected.re <- matrix(nrow = 4, ncol = n.scenarios, dimnames = disconnected.summary.matrix.names )
  coverage.summary.disconnected.re <- coverage.mean.disconnected.re <- coverage.ll.disconnected.re <- coverage.ul.disconnected.re  <- matrix(nrow = 4, ncol = n.scenarios, dimnames = disconnected.summary.matrix.names )
  
  for(i.scenario in 1:n.scenarios)
  {
    # Mean and 95% limits for use in figure
    bias.mean.connected.re[, i.scenario] <- colMeans(bias.tables.re[[i.scenario]], na.rm = TRUE)[1:5]
    bias.mean.disconnected.re[, i.scenario] <- colMeans(bias.tables.re[[i.scenario]], na.rm = TRUE)[6:9]
    coverage.mean.connected.re[, i.scenario] <- colMeans(coverage.tables.re[[i.scenario]], na.rm = TRUE)[1:5]
    coverage.mean.disconnected.re[, i.scenario] <- colMeans(coverage.tables.re[[i.scenario]], na.rm = TRUE)[6:9]
    
    bias.ll.connected.re[, i.scenario] <- apply(bias.tables.re[[i.scenario]], c(2), mean)[1:5] - 1.96 * apply(bias.tables.re[[i.scenario]], c(2), sd)[1:5]/sqrt(n.simulations)
    bias.ll.disconnected.re[, i.scenario] <- apply(bias.tables.re[[i.scenario]], c(2), mean)[6:9] - 1.96 * apply(bias.tables.re[[i.scenario]], c(2), sd)[6:9]/sqrt(n.simulations)
    coverage.ll.connected.re[, i.scenario] <-  apply(coverage.tables.re[[i.scenario]], c(2), mean)[1:5] - 1.96 * apply(coverage.tables.re[[i.scenario]], c(2), sd)[1:5]/sqrt(n.simulations)
    coverage.ll.disconnected.re[, i.scenario] <-  apply(coverage.tables.re[[i.scenario]], c(2), mean)[6:9] - 1.96 * apply(coverage.tables.re[[i.scenario]], c(2), sd)[6:9]/sqrt(n.simulations)
    
    bias.ul.connected.re[, i.scenario] <- apply(bias.tables.re[[i.scenario]], c(2), mean)[1:5] + 1.96 * apply(bias.tables.re[[i.scenario]], c(2), sd)[1:5]/sqrt(n.simulations)
    bias.ul.disconnected.re[, i.scenario] <- apply(bias.tables.re[[i.scenario]], c(2), mean)[6:9] + 1.96 * apply(bias.tables.re[[i.scenario]], c(2), sd)[6:9]/sqrt(n.simulations)
    coverage.ul.connected.re[, i.scenario] <-  apply(coverage.tables.re[[i.scenario]], c(2), mean)[1:5] + 1.96 * apply(coverage.tables.re[[i.scenario]], c(2), sd)[1:5]/sqrt(n.simulations)
    coverage.ul.disconnected.re[, i.scenario] <-  apply(coverage.tables.re[[i.scenario]], c(2), mean)[6:9] + 1.96 * apply(coverage.tables.re[[i.scenario]], c(2), sd)[6:9]/sqrt(n.simulations)
    
    
  }
  # Ensure coverage does not extend beyond 0, 1
  coverage.ll.connected.re[coverage.ll.connected.re < 0] <- 0
  coverage.ll.disconnected.re[coverage.ll.disconnected.re < 0] <- 0
  coverage.ul.connected.re[coverage.ul.connected.re > 1] <- 1
  coverage.ul.disconnected.re[coverage.ul.disconnected.re > 1] <- 1
  for(i.scenario in 1:n.scenarios) 
  {
    # Summary for use in tables
    bias.summary.connected.re[, i.scenario] <- paste0(format(bias.mean.connected.re[, i.scenario], digits = 2), " (",
                                                      format(bias.ll.connected.re[, i.scenario], digits = 2), ", ",
                                                      format(bias.ul.connected.re[, i.scenario], digits = 2), ")")
    
    bias.summary.disconnected.re[, i.scenario] <- paste0(format(bias.mean.disconnected.re[, i.scenario], digits = 2), " (",
                                                         format(bias.ll.disconnected.re[, i.scenario], digits = 2), ", ",
                                                         format(bias.ul.disconnected.re[, i.scenario], digits = 2), ")")
    coverage.summary.connected.re[, i.scenario] <- paste0(format(coverage.mean.connected.re[, i.scenario], digits = 2), " (",
                                                          format(coverage.ll.connected.re[, i.scenario], digits = 2), ", ",
                                                          format(coverage.ul.connected.re[, i.scenario], digits = 2), ")")
    coverage.summary.disconnected.re[, i.scenario] <- paste0(format(coverage.mean.disconnected.re[, i.scenario], digits = 2), " (",
                                                             format(coverage.ll.disconnected.re[, i.scenario], digits = 2), ", ",
                                                             format(coverage.ul.disconnected.re[, i.scenario], digits = 2), ")")
    
  }
  
  plot.simulation.results(bias.mean.connected.re,bias.mean.disconnected.re,bias.ll.connected.re, bias.ll.disconnected.re, bias.ul.connected.re, bias.ul.disconnected.re,
                                 coverage.mean.connected.re,coverage.mean.disconnected.re,coverage.ll.connected.re, coverage.ll.disconnected.re, coverage.ul.connected.re, coverage.ul.disconnected.re,
                          start.of.filename = paste0(baseline.directory,"/results/plots/simulation.re.",n.simulations, ".", n.patients, "patients"))
  
                          
  
  save(bias.summary.connected.re,bias.summary.disconnected.re,coverage.summary.connected.re,coverage.summary.disconnected.re,bias.tables.re,coverage.tables.re,
       bias.mean.connected.re,bias.mean.disconnected.re,bias.ll.connected.re, bias.ll.disconnected.re, bias.ul.connected.re, bias.ul.disconnected.re,
       coverage.mean.connected.re,coverage.mean.disconnected.re,coverage.ll.connected.re, coverage.ll.disconnected.re, coverage.ul.connected.re, coverage.ul.disconnected.re,
       file=paste0(baseline.directory,"/code/simulation.results.re.",n.simulations,".", n.patients, "patients.v3.rda"))
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "Bias conn")
  writeData(wb, "Bias conn", cbind(rownames(bias.summary.connected.re), bias.summary.connected.re), startRow = 1, startCol = 1)
  addWorksheet(wb, "Bias disc")
  writeData(wb, "Bias disc", cbind(rownames(bias.summary.disconnected.re), bias.summary.disconnected.re), startRow = 1, startCol = 1)
  addWorksheet(wb, "Coverage conn")
  writeData(wb, "Coverage conn", cbind(rownames(coverage.summary.connected.re), coverage.summary.connected.re), startRow = 1, startCol = 1)
  addWorksheet(wb, "Coverage disc")
  writeData(wb, "Coverage disc", coverage.summary.disconnected.re, startRow = 1, startCol = 2)
  saveWorkbook(wb, file = paste0(baseline.directory,"/results/simulation.results.re.",n.simulations,".", n.patients, "patients.v3.xlsx"), overwrite = TRUE)    
  
  
  
}




