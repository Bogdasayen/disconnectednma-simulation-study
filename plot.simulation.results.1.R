plot.simulation.results <- function(bias.mean.connected,bias.mean.disconnected,bias.ll.connected, bias.ll.disconnected, bias.ul.connected, bias.ul.disconnected,
                        coverage.mean.connected,coverage.mean.disconnected,coverage.ll.connected, coverage.ll.disconnected, coverage.ul.connected, coverage.ul.disconnected,
                        start.of.filename)
{
  n.methods = dim(bias.mean.connected)[1]
  n.scenarios = dim(bias.mean.connected)[2]
  par(mfrow = c(1,1))
  
  # Plot bias in connected network
  win.metafile(paste0(start.of.filename, ".bias.connected.wmf"))
  plot(bias.mean.connected[1,], ylim = c(min(bias.ll.connected), max(bias.ul.connected)), xaxt = "n", xlab = "Scenario", ylab = "")
  title("Bias connected")
  axis(side = 1, at = c(1:n.scenarios), labels = c(1:12))
  legend("topright", col = c(1:n.methods), lty = 1, legend = c("RCT only", "ALM single-arms", "ALM disconnected RCTs", "RP single-arms", "RP disconnected RCTs"))
  for(i.scenario in 1:n.scenarios) {
    lines(x= c(i.scenario, i.scenario), y = c(bias.ll.connected[1, i.scenario], bias.ul.connected[1, i.scenario]))
  }
  for(i.method in 2:n.methods) {
    points(y = bias.mean.connected[i.method, ], x = c(1:n.scenarios) + i.method/10, col = i.method)
    for(i.scenario in 1:n.scenarios) {
      lines(x= c(i.scenario, i.scenario) + i.method/10, y = c(bias.ll.connected[i.method, i.scenario], bias.ul.connected[i.method, i.scenario]), col = i.method)
    }
  }
  dev.off()
  
  # Plot bias in disconnected network
  win.metafile(paste0(start.of.filename, ".bias.disconnected.wmf"))
  plot(bias.mean.disconnected[1,], ylim = c(min(bias.ll.disconnected), max(bias.ul.disconnected)), xaxt = "n", xlab = "Scenario", ylab = "", col = 2)
  title("Bias disconnected")
  axis(side = 1, at = c(1:n.scenarios), labels = c(1:12))
  for(i.scenario in 1:n.scenarios) {
    lines(x= c(i.scenario, i.scenario), y = c(bias.ll.disconnected[1, i.scenario], bias.ul.disconnected[1, i.scenario]), col = 2)
  }
  for(i.method in 2:(n.methods-1)) {
    points(y = bias.mean.disconnected[i.method, ], x = c(1:n.scenarios) + i.method/10, col = i.method + 1)
    for(i.scenario in 1:n.scenarios) {
      lines(x= c(i.scenario, i.scenario) + i.method/10, y = c(bias.ll.disconnected[i.method, i.scenario], bias.ul.disconnected[i.method, i.scenario]), col = i.method + 1)
    }
  }
  dev.off()
  
  # Plot coverage in connected network
  win.metafile(paste0(start.of.filename, ".coverage.connected.wmf"))
  plot(coverage.mean.connected[1,], ylim = c(min(coverage.ll.connected), max(coverage.ul.connected)), xaxt = "n", xlab = "Scenario", ylab = "")
  title("Coverage connected")
  axis(side = 1, at = c(1:n.scenarios), labels = c(1:12))
  for(i.scenario in 1:n.scenarios) {
    lines(x= c(i.scenario, i.scenario), y = c(coverage.ll.connected[1, i.scenario], coverage.ul.connected[1, i.scenario]))
  }
  for(i.method in 2:n.methods) {
    points(y = coverage.mean.connected[i.method, ], x = c(1:n.scenarios) + i.method/10, col = i.method)
    for(i.scenario in 1:n.scenarios) {
      lines(x= c(i.scenario, i.scenario) + i.method/10, y = c(coverage.ll.connected[i.method, i.scenario], coverage.ul.connected[i.method, i.scenario]), col = i.method)
    }
  }
  dev.off()
  
  # Plot coverage in disconnected network
  win.metafile(paste0(start.of.filename, ".coverage.disconnected.wmf"))
  plot(coverage.mean.disconnected[1,], ylim = c(min(coverage.ll.disconnected), max(coverage.ul.disconnected)), xaxt = "n", xlab = "Scenario", ylab = "", col = 2)
  title("Coverage disconnected")
  axis(side = 1, at = c(1:n.scenarios), labels = c(1:12))
  for(i.scenario in 1:n.scenarios) {
    lines(x= c(i.scenario, i.scenario), y = c(coverage.ll.disconnected[1, i.scenario], coverage.ul.disconnected[1, i.scenario]), col = 2)
  }
  for(i.method in 2:(n.methods-1)) {
    points(y = coverage.mean.disconnected[i.method, ], x = c(1:n.scenarios) + i.method/10, col = i.method + 1)
    for(i.scenario in 1:n.scenarios) {
      lines(x= c(i.scenario, i.scenario) + i.method/10, y = c(coverage.ll.disconnected[i.method, i.scenario], coverage.ul.disconnected[i.method, i.scenario]), col = i.method + 1)
    }
  }
  dev.off()
  
}