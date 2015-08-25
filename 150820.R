# install.packages("adaptivetau")
# install.packages("ggplot2")
# install.packages("RColorBrewer")
# install.packages("lattice")
# install.packages("deSolve")
# install.packages("reshape2")
# install.packages("gplots")
library(adaptivetau)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(deSolve)
library(reshape2)
library(gplots)

# hello gitword #

######################
# Plotting Functions #
######################
# Compartments #
Plot_Compartments <- function(trial){  
  ggplot(data=trial, aes(x=time)) +
    geom_line(aes(x=time, y=S_0), col="blue", linetype="dotted") +
    geom_line(aes(x=time, y=S_1), col="blue", linetype="dashed") +
    geom_line(aes(x=time, y=S_2), col="blue", linetype="solid") +
    geom_line(aes(x=time, y=I_0), col="red", linetype="dotted") +
    geom_line(aes(x=time, y=I_1), col="red", linetype="dashed") +
    geom_line(aes(x=time, y=I_2), col="red", linetype="solid") +
    geom_line(aes(x=time, y=R_0), col="green", linetype="dotted") +
    geom_line(aes(x=time, y=R_1), col="green", linetype="dashed") +
    geom_line(aes(x=time, y=R_2), col="green", linetype="solid") +
    theme_bw() + ylab("Count") + ggtitle("Compartments") 
}
# Plot_Compartments(trial)
# Plot_Compartments(simulation.df)

# SEIR compiled #
Plot_SEIR <- function(trial){
  ggplot(data=trial, aes(x=time)) +
    geom_line(aes(x=time, y=(S_0 + S_1 + S_2)), col="blue") +
    geom_line(aes(x=time, y=(I_0 + I_1 + I_2)), col="red") +
    geom_line(aes(x=time, y=(R_0 + R_1 + R_2)), col="green") +
    theme_bw() + ylab("Count") + ggtitle("SEIR")
}
# Plot_SEIR(trial)
# Plot_SEIR(simulation.df)

####################
# Stochastic Model #
####################
## x is a vector(S_0, S_1, S_2, E_0, E_1, E_2, I_0, I_1, I_2, R_0, R_1, R_2)
SEIRrates<-function(x,params,t){
  with(params,{
    return(c(b*ntot,    ## birth
             u*x[1],    ## death of S_0
             u*x[2],    ## death of S_1
             u*x[3],    ## death of S_2
             u*x[4],    ## death of E_0
             u*x[5],    ## death of E_1
             u*x[6],    ## death of E_2
             u*x[7],    ## death of I_0
             u*x[8],    ## death of I_1
             u*x[9],    ## death of I_2
             u*x[10],  	## death of R_0
             u*x[11],  	## death of R_1
             u*x[12],  	## death of R_2
             (B0*x[1]*x[7] + (1-k1)*B0*x[1]*x[8] + (1-k2)*B0*x[1]*x[9])/ntot,           ## infection of unvaccinated susceptible person
             (1-vac1)*(B0*x[2]*x[7] + (1-k1)*B0*x[2]*x[8] + (1-k2)*B0*x[2]*x[9])/ntot,  ## infection of once-vaccinated susceptible person
             (1-vac2)*(B0*x[3]*x[7] + (1-k1)*B0*x[3]*x[8] + (1-k2)*B0*x[3]*x[9])/ntot,  ## infection of twice-vaccinated susceptible person
             e*x[4],     ## I_0 becomes infectious
               e*x[5],     ## I_1 becomes infectious
             e*x[6],     ## I_2 becomes infectious
             r*x[7],    ## recovery of I_0
             r*x[8],    ## recovery of I_1
             r*x[9]     ## recovery of I_2
    ))
  })
}

#############################################
# Paramaters and inits for Stochastic Model #
#############################################
par=list(u = 0,            # Death rate
         b = 0,            # Birth rate
         B0 = 0.5,         # Transmission parameter for a non-vaccinated susceptible person contacting a non-vaccinated infectious person
         k1 = 0.2,         # Proportional infectiousness of a once-vaccinated person
         k2 = 0.5,         # Proportional infectiousness of a twice-vaccinated person
         vac1 = 0.5,       # Personal efficacy of one dose of vaccine
         vac2 = 0.8,       # Personal efficacy of two doses of vaccine
         r = 1/3,          # Recovery rate (1/days)
         e = 1/2,          # Incubation period
         ntot = 25000,     # Population size
         v1_day = 30,      # Day of first vaccine dose (If you don't want vaccination, then set this to 100 or whatever your t_final is)
         v_delay = 14,     # Delay until second vaccine dose
         v_count = 2000)   # Total number of vaccines available
inits <- c(par$ntot-10,0,0,0,0,0,10,0,0,0,0,0)
inits_original <- inits
par_original <- par # Save original pars for later use
campaign <- 0

#########################
# Run stochastic models #
#########################
Run_Model_S <- function (par, inits, tf=365, campaign){
  # determine distribution of vaccines
  if (campaign == 0){
    v1_count = 0
    v2_count = 0
  } else if (campaign == 1){
    v1_count = par$v_count
    v2_count = 0
  } else if (campaign == 2){
    v1_count = par$v_count/2
    v2_count = par$v_count/2
  }
  # section time depending on campaign
  segment_0 = data.frame(                ## segment_0 is before the first dose campaign on v1_day
    ssa.adaptivetau(inits,
                    matrix(c(
                      1,0,0,0,0,0,0,0,0,0,0,0,   ## birth
                      -1,0,0,0,0,0,0,0,0,0,0,0,  ## death...
                      0,-1,0,0,0,0,0,0,0,0,0,0,
                      0,0,-1,0,0,0,0,0,0,0,0,0,
                      0,0,0,-1,0,0,0,0,0,0,0,0,
                      0,0,0,0,-1,0,0,0,0,0,0,0,
                      0,0,0,0,0,-1,0,0,0,0,0,0,
                      0,0,0,0,0,0,-1,0,0,0,0,0,
                      0,0,0,0,0,0,0,-1,0,0,0,0,
                      0,0,0,0,0,0,0,0,-1,0,0,0,
                      0,0,0,0,0,0,0,0,0,-1,0,0,
                      0,0,0,0,0,0,0,0,0,0,-1,0,
                      0,0,0,0,0,0,0,0,0,0,0,-1,  ## ...death
                      -1,0,0,1,0,0,0,0,0,0,0,0,  ## infection...
                      0,-1,0,0,1,0,0,0,0,0,0,0,
                      0,0,-1,0,0,1,0,0,0,0,0,0,  ## ...infection
                      0,0,0,-1,0,0,1,0,0,0,0,0,  ## infectiousness
                      0,0,0,0,-1,0,0,1,0,0,0,0,
                      0,0,0,0,0,-1,0,0,1,0,0,0,  ## ...infectiousness
                      0,0,0,0,0,0,-1,0,0,1,0,0,  ## recovery...
                      0,0,0,0,0,0,0,-1,0,0,1,0,
                      0,0,0,0,0,0,0,0,-1,0,0,1   ## ...recovery
                    ),nrow=12),
                    SEIRrates,
                    par,
                    tf=par$v1_day,
                    tl.params=list(epsilon=0.02)) )
  names(segment_0) <- c("time", "S_0", "S_1", "S_2", "E_0", "E_1", "E_2", "I_0", "I_1", "I_2", "R_0", "R_1", "R_2")
  
  if (sum(segment_0[nrow(segment_0), c("E_0","E_1","E_2","I_0", "I_1", "I_2")]) > 0 &             ## If you still have infected people
        segment_0[nrow(segment_0),"time"] < tf){                               ## If you haven't reached the 100th day
    if (segment_0[nrow(segment_0),"S_0"] > v1_count){                       ## If you enough S to vaccinate v1_count
      count = v1_count                                                      ## then you'll vaccinate v_count
    } else {count = segment_0[nrow(segment_0),"S_0"]}                           ## otherwise, you'll vaccinate all the S_0 remaining
    
    inits_1 <- as.numeric(segment_0[nrow(segment_0),c("S_0", "S_1", "S_2", "E_0", "E_1", "E_2", "I_0", "I_1", "I_2", "R_0", "R_1", "R_2")])
    inits_1[1] <- as.numeric(inits_1[1] - count)
    inits_1[2] <- as.numeric(inits_1[2] + count)
    
    segment_1 = data.frame(                  ## segment_1 is between the first dose campaign on v1_day and the second dose v_delay later
      ssa.adaptivetau(inits_1,
                      matrix(c(
                        1,0,0,0,0,0,0,0,0,0,0,0,   ## birth
                        -1,0,0,0,0,0,0,0,0,0,0,0,  ## death...
                        0,-1,0,0,0,0,0,0,0,0,0,0,
                        0,0,-1,0,0,0,0,0,0,0,0,0,
                        0,0,0,-1,0,0,0,0,0,0,0,0,
                        0,0,0,0,-1,0,0,0,0,0,0,0,
                        0,0,0,0,0,-1,0,0,0,0,0,0,
                        0,0,0,0,0,0,-1,0,0,0,0,0,
                        0,0,0,0,0,0,0,-1,0,0,0,0,
                        0,0,0,0,0,0,0,0,-1,0,0,0,
                        0,0,0,0,0,0,0,0,0,-1,0,0,
                        0,0,0,0,0,0,0,0,0,0,-1,0,
                        0,0,0,0,0,0,0,0,0,0,0,-1,  ## ...death
                        -1,0,0,1,0,0,0,0,0,0,0,0,  ## infection...
                        0,-1,0,0,1,0,0,0,0,0,0,0,
                        0,0,-1,0,0,1,0,0,0,0,0,0,  ## ...infection
                        0,0,0,-1,0,0,1,0,0,0,0,0,  ## infectiousness
                        0,0,0,0,-1,0,0,1,0,0,0,0,
                        0,0,0,0,0,-1,0,0,1,0,0,0,  ## ...infectiousness
                        0,0,0,0,0,0,-1,0,0,1,0,0,  ## recovery...
                        0,0,0,0,0,0,0,-1,0,0,1,0,
                        0,0,0,0,0,0,0,0,-1,0,0,1   ## ...recovery
                      ),nrow=12),
                      SEIRrates,
                      par,
                      tf=par$v_delay,       ## Note that the timer starts at zero, so we only go v_delay more days in segment_1
                      tl.params=list(epsilon=0.02)) )
    names(segment_1) <- c("time", "S_0", "S_1", "S_2", "E_0", "E_1", "E_2", "I_0", "I_1", "I_2", "R_0", "R_1", "R_2")
    segment_1$time <- segment_1$time + par$v1_day
    
    if (sum(segment_1[nrow(segment_1), c("E_0","E_1","E_2","I_0", "I_1", "I_2")]) > 0 &            ## If you still have infected people
          segment_1[nrow(segment_1),"time"] < tf){                                ## If you haven't reached the 100th day)  
      if (segment_1[nrow(segment_1),"S_1"] > v2_count){                      ## If you enough S_1 people to vaccinate v2_count of them
        count = v2_count                                                     ## then vaccinate v2_count
        count_extra = 0                                                          ## and you have no extra vaccines to give to S_0
      } else {
        count = segment_1[nrow(segment_1),"S_1"]                                 ## If you do have leftover vaccines
        count_extra = v2_count - count                                       
      } 
      
      inits_2 <- as.numeric(segment_1[nrow(segment_1),c("S_0", "S_1", "S_2", "E_0", "E_1", "E_2", "I_0", "I_1", "I_2", "R_0", "R_1", "R_2")])
      inits_2[2] <- as.numeric(inits_2[2] - count)
      inits_2[3] <- as.numeric(inits_2[3] + count)
      inits_2[2] <- as.numeric(inits_2[2] + inits_2[1] - max(0, as.numeric(inits_2[1] - count_extra)) )    ## any extra vaccines will bring up to count_extra people from S_0 to S_1
      inits_2[1] <- max(0, as.numeric(inits_2[1] - count_extra))
      
      segment_2 = data.frame(                ## segment_2 is after the second dose campaign on v1_day+v_delay
        ssa.adaptivetau(inits_2,
                        matrix(c(
                          1,0,0,0,0,0,0,0,0,0,0,0,   ## birth
                          -1,0,0,0,0,0,0,0,0,0,0,0,  ## death...
                          0,-1,0,0,0,0,0,0,0,0,0,0,
                          0,0,-1,0,0,0,0,0,0,0,0,0,
                          0,0,0,-1,0,0,0,0,0,0,0,0,
                          0,0,0,0,-1,0,0,0,0,0,0,0,
                          0,0,0,0,0,-1,0,0,0,0,0,0,
                          0,0,0,0,0,0,-1,0,0,0,0,0,
                          0,0,0,0,0,0,0,-1,0,0,0,0,
                          0,0,0,0,0,0,0,0,-1,0,0,0,
                          0,0,0,0,0,0,0,0,0,-1,0,0,
                          0,0,0,0,0,0,0,0,0,0,-1,0,
                          0,0,0,0,0,0,0,0,0,0,0,-1,  ## ...death
                          -1,0,0,1,0,0,0,0,0,0,0,0,  ## infection...
                          0,-1,0,0,1,0,0,0,0,0,0,0,
                          0,0,-1,0,0,1,0,0,0,0,0,0,  ## ...infection
                          0,0,0,-1,0,0,1,0,0,0,0,0,  ## infectiousness
                          0,0,0,0,-1,0,0,1,0,0,0,0,
                          0,0,0,0,0,-1,0,0,1,0,0,0,  ## ...infectiousness
                          0,0,0,0,0,0,-1,0,0,1,0,0,  ## recovery...
                          0,0,0,0,0,0,0,-1,0,0,1,0,
                          0,0,0,0,0,0,0,0,-1,0,0,1   ## ...recovery
                        ),nrow=12),
                        SEIRrates,
                        par,
                        tf=tf - par$v1_day - par$v_delay,
                        tl.params=list(epsilon=0.02)) )
      names(segment_2) <- c("time", "S_0", "S_1", "S_2", "E_0", "E_1", "E_2", "I_0", "I_1", "I_2", "R_0", "R_1", "R_2")
      segment_2$time <- segment_2$time + par$v1_day + par$v_delay
      trial = rbind(segment_0, segment_1) 
      trial = rbind(trial, segment_2)
      
    } else{trial = rbind(segment_0, segment_1)}
  } else{trial = segment_0}
  
  trial <- data.frame(trial)
  names(trial) <- c("time", "S_0", "S_1", "S_2", "E_0", "E_1", "E_2", "I_0", "I_1", "I_2", "R_0", "R_1", "R_2")
  return(trial)
}

#########################
# Stochastic Simulation #
#########################
trial <- Run_Model_S(par=par, inits=inits, campaign=1)
Plot_Compartments(trial)
Plot_SEIR(trial)

#######################
# Deterministic Model #
#######################
Model_D <- function(t, x, parms){ 
  with(as.list(c(parms,x)),{
    dS_0    <- -(B0*S_0*I_0/N + (1-k1)*B0*S_0*I_1/N + (1-k2)*B0*S_0*I_2/N)
    dS_1  <- -(1-vac1)*(B0*S_1*I_0/N + (1-k1)*B0*S_1*I_1/N + (1-k2)*B0*S_1*I_2/N)
    dS_2  <- -(1-vac2)*(B0*S_2*I_0/N + (1-k1)*B0*S_2*I_1/N + (1-k2)*B0*S_2*I_2/N)
    dE_0  <- (B0*S_0*I_0/N + (1-k1)*B0*S_0*I_1/N + (1-k2)*B0*S_0*I_2/N) - e*E_0
    dE_1  <- (1-vac1)*(B0*S_1*I_0/N + (1-k1)*B0*S_1*I_1/N + (1-k2)*B0*S_1*I_2/N) - e*E_1
    dE_2  <- (1-vac2)*(B0*S_2*I_0/N + (1-k1)*B0*S_2*I_1/N + (1-k2)*B0*S_2*I_2/N) - e*E_2
    dI_0  <- e*E_0 - r*I_0
    dI_1  <- e*E_1 - r*I_1
    dI_2  <- e*E_2 - r*I_2
    dR_0  <- r*I_0
    dR_1  <- r*I_1
    dR_2  <- r*I_2
    dN    <- 0
    der   <- c(dS_0, dS_1, dS_2, dE_0, dE_1, dE_2, dI_0, dI_1, dI_2, dR_0, dR_1, dR_2, dN)
    list(der) #output
  })
}

##################################################
# Vaccination conditions for deterministic model #
##################################################
Run_Model_D <- function(inits, dt, parms){  
########No vaccine
  if (parms["v1_count"] == 0 && parms["v2_count"] == 0){
    simulation <- data.frame(lsoda(inits, dt, Model_D, parms=parms))
  } else if (parms["v2_count"] == 0){   
########Vaccinate with one dose only
    dt_0 <- dt[dt < parms["v1_day"]]
    dt_1 <- dt[dt >= parms["v1_day"]]
    segment_0 <- data.frame(lsoda(inits, dt_0, Model_D, parms=parms))
    segment_1_inits <- unlist(segment_0[nrow(segment_0),c("S_0", "S_1", "S_2", "E_0", "E_1", "E_2", "I_0", "I_1", "I_2", "R_0", "R_1", "R_2","N")])
    #Update inits by vaccinating
    if (segment_1_inits[1] < parms["v1_count"]){  #If you don't have enough S_0 to vaccinate
      segment_1_inits[2] <- segment_1_inits[1]
      segment_1_inits[1] <- 0
    } else {
      segment_1_inits[2] <- parms["v1_count"]
      segment_1_inits[1] <- segment_1_inits[1] - parms["v1_count"]
    }
    segment_1 <- data.frame(lsoda(segment_1_inits,dt_1, Model_D, parms=parms))
    simulation <- rbind(segment_0, segment_1)
  } else {
########Vaccinate with two doses
    dt_0 <- dt[dt <  parms["v1_day"]]
    dt_1 <- dt[dt >= parms["v1_day"] & dt < parms["v2_day"]]
    dt_2 <- dt[dt >= parms["v2_day"]]
    segment_0 <- data.frame(lsoda(inits, dt_0, Model_D, parms=parms))
    segment_1_inits <- unlist(segment_0[nrow(segment_0),c("S_0", "S_1", "S_2", "E_0", "E_1", "E_2", "I_0", "I_1", "I_2", "R_0", "R_1", "R_2","N")])
    #Update inits by vaccinating
    if (segment_1_inits[1] < parms["v1_count"]){  #If you don't have enough S_1 to vaccinate
      segment_1_inits[2] <- segment_1_inits[1]
      segment_1_inits[1] <- 0
    } else {
      segment_1_inits[2] <- parms["v1_count"]
      segment_1_inits[1] <- segment_1_inits[1] - parms["v1_count"]
    }
    segment_1 <- data.frame(lsoda(segment_1_inits,dt_1, Model_D, parms=parms))
    segment_2_inits <- unlist(segment_1[nrow(segment_1),c("S_0", "S_1", "S_2", "E_0", "E_1", "E_2", "I_0", "I_1", "I_2", "R_0", "R_1", "R_2","N")])
    #Update inits by vaccinating again
    if (segment_2_inits[2] < parms["v2_count"]){  #If you don't have enough S_1 to vaccinate
      S1_to_S2 = segment_2_inits[2]                               ## If you do have leftover vaccines
      count_extra = parms["v2_count"] - S1_to_S2
      S0_to_S1 = min(count_extra, segment_2_inits[1])
      segment_2_inits[3] <- S1_to_S2
      segment_2_inits[2] <- segment_2_inits[2] - S1_to_S2 + S0_to_S1
      segment_2_inits[1] <- segment_2_inits[1] - S0_to_S1
    } else {     # You do have enough S_1 to vaccinate again
      segment_2_inits[3] <- parms["v2_count"]
      segment_2_inits[2] <- segment_2_inits[2] - parms["v2_count"]
    }      
    segment_2 <- data.frame(lsoda(segment_2_inits,dt_2, Model_D, parms=parms))
    simulation <- rbind(segment_0, segment_1)
    simulation <- rbind(simulation, segment_2)
  }
}

################################################
# Paramaters and inits for Deterministic Model #
################################################
parms = c(u = 0,            # Death rate
          b = 0,            # Birth rate
          B0 = 0.5,         # Transmission parameter for a non-vaccinated susceptible person contacting a non-vaccinated infectious person
          k1 = 0.0,         # Reduction in infectiousness of a once-vaccinated person
          k2 = 0.0,         # Reduction in infectiousness of a twice-vaccinated person
          vac1 = 0.5,       # Personal efficacy of one dose of vaccine
          vac2 = 0.85,      # Personal efficacy of two doses of vaccine
          r = 1/3,          # Recovery rate (1/days)
          e = 1/2,          # Incubation (1/days)
          ntot = 25000,     # Population size
          v1_day = 50,      # Day of first vaccine dose (If you don't want vaccination, then set this to 100 or whatever your t_final is)
          v2_day = 64,      # Day of second vaccine dose
          v1_count = 1000,  # number of vaccines intended for first dose
          v2_count = 1000)  # number of vaccines intended for second dose 

# save parms for later use
parms_original <- parms

t_final <- 365
dt <-seq(from=1, to=t_final, by=1)

inits <- c(S_0=25000-1,
           S_1=0,
           S_2=0,
           E_0=0,
           E_1=0,
           E_2=0,
           I_0=1,
           I_1=0,
           I_2=0,
           R_0=0,
           R_1=0,
           R_2=0,
           N  =25000)
inits_original <- inits

############################
# Deterministic Simulation #
############################
parms["ntot"] <- 25000
parms["B0"] <- .60
simulation <- Run_Model_D(inits, dt, parms=parms)
simulation.df <- as.data.frame(simulation)
simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
I_sum <- simulation.df[c("I_0")] + simulation.df[c("I_1")] + simulation.df[c("I_2")]
peak <- which.max(I_sum[,1])
#Plot_Compartments(simulation.df)
Plot_SEIR(simulation.df)
Rnaught <- parms["B0"]/parms["r"]
CAR <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])/parms["ntot"]
Rnaught
CAR

#######################
# Heat Map Comparison #
#######################
bins <- 10
repeats <- 50

# create a red color palette 
jOrRdFun <- colorRampPalette(brewer.pal(n = 9, "OrRd"))
paletteSize <- 256
jOrRdPalette <- jOrRdFun(paletteSize)

# blue color palette
jBluesFun <- colorRampPalette(brewer.pal(n = 9, "Blues"))
jBluesPalette <- jBluesFun(paletteSize)

# green color palette
jGreensFun <- colorRampPalette(brewer.pal(n = 9, "Greens"))
jGreensPalette <- jGreensFun(paletteSize)

# red to blue color palette
jRdBuFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))
jRdBuPalette <- jRdBuFun(paletteSize)

# red to green color palette
jRdYlGnFun <- colorRampPalette(brewer.pal(n = 9, "RdYlGn"))
jRdYlGnPalette <- jRdYlGnFun(paletteSize)

# blue to green color palette
jPRGnFun <- colorRampPalette(brewer.pal(n = 9, "PRGn"))
jPRGnPalette <- jPRGnFun(paletteSize)
jBluGnPalette <- jPRGnPalette
for (i in 1:128){
  jBluGnPalette[i] <- jRdBuPalette[257-i]
}
jBluGrnPalette <- jBluGnPalette
for (i in 1:256){
  jBluGrnPalette[i] <- jBluGnPalette[257-i]
}

##################################################
##################################################
## end required parameters, functions, palettes ##
##################################################
##################################################

############################
# Experiment 1: Vary inits #
############################
# No vax
# trial <- Run_Model_S(par=par, inits = c(990,0,0,0,0,0,10,0,0,0,0,0))
# print(sum(trial[nrow(trial), c("R_0", "R_1","R_2")]))

# # All once-vax
# trial <- Run_Model_S(par=par, inits = c(0,990,0,0,0,0,0,10,0,0,0,0))
# print(sum(trial[nrow(trial), c("R_0", "R_1","R_2")]))
# 
# # All twice-vax
# trial <- Run_Model_S(par=par, inits = c(0,0,990,0,0,0,0,0,10,0,0,0))
# print(sum(trial[nrow(trial), c("R_0", "R_1","R_2")]))

# Even mix. Check for dose response
par$B0 = .25
trial <- Run_Model_S(par=par, inits = c(9999,0,0,0,0,0,1,0,0,0,0,0), campaign=2)
par$v1_day = 10
Plot_Compartments(trial)

# Efficacious vaccine after 15 days, 10 day second dose delay
par$B0 = .3
par$v1_day = 15;
par$v_delay = 10;
trial <- Run_Model_S(par=par, inits = c(9990,0,0,0,0,0,10,0,0,0,0,0), tf=300, campaign=2)
print(sum(trial[nrow(trial), c("R_0", "R_1","R_2")]))
Plot_Compartments(trial)

# normal epidemic
parms <- parms_original
parms["B0"] <- .33
parms["r"] <- 1/4
simulation <- Run_Model_D(inits, dt, parms=parms)
simulation.df <- as.data.frame(simulation)
simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
Plot_SEIR(simulation.df)

#######################################
# Experiment 2: demonstration of loop #
#######################################
reps=100
Exp2 <- data.frame(1:length(reps))
Exp2[,2] <- NA
names(Exp2) <- c("Trial","Count")
for (i in 1:reps){
  trial <- Run_Model_S(par=par, inits = c(9990,0,0,0,0,0,10,0,0,0,0,0))
  Exp2[i,1] <- i
  Exp2[i,2] <- sum(trial[nrow(trial), c("R_0", "R_1","R_2")])    # Count the number of people who recovered from the disease by day 100
  print(i)
}
View(Exp2)

######################################################
# Experiment 3: Compare Stochastic and Deterministic #
######################################################

comparisons <- 10
compare <- data.frame(c(rep(0,comparisons)))
for (i in 1:comparisons){
  inits_s <- c(9999,0,0,0,0,0,1,0,0,0,0,0)
  trial_s <- Run_Model_S(par = par, inits = inits_s, campaign = 0)
  compare[i,1] <- sum(trial_s[nrow(trial_s), c("R_0", "R_1", "R_2")])
}
hist(compare[,1])
summary(compare[compare[,1]>1000,1])
avg_s <- sum(compare)/nrow(compare)
inits_d <- c(S_0=0.9999,
             S_1=0.00,
             S_2=0.00,
             E_0=0.00,
             E_1=0.00,
             E_2=0.00,
             I_0=0.0001,
             I_1=0.00,
             I_2=0.00,
             R_0=0.00,
             R_1=0.00,
             R_2=0.00)
simulation <- lsoda(inits_d, dt, Run_Model_D, parms=parms)
simulation.df <- as.data.frame(simulation)
simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
trial_d <- 10000*simulation.df
avg_d <- sum(trial_d[nrow(trial_d), c("R_0", "R_1", "R_2")])

###############################################
# Compare infectiousness of vaccinated people #
###############################################
par <- par_original # reset parameters
par$vac1 = 0.5   # set the protective efficacy of one or two doses both to 0.5
par$vac2 = 0.5 

heat_map_k <- data.frame(matrix(rep(0, bins*bins), nrow=bins))

for (x in 1:bins){             # changes in k1
  par$k1 = .1*x
  rownames(heat_map_k)[x] <- as.character(par$k1)
  for (y in 1:x){              # changes in k2
    par$k2 = .1*y
    colnames(heat_map_k)[y] <- as.character(par$k2)
    for (z in 1:repeats){      # run trials and average number recovered
      trial <- Run_Model_S(par=par, inits = c(3300,3300,3300,0,0,0,40,30,30,0,0,0), campaign = 0)
      heat_map_k[x,y] <- heat_map_k[x,y] + sum(trial[nrow(trial), c("R_0", "R_1","R_2")])
    }
  }
}
heat_map_k <- heat_map_k/repeats

# create a heatmap 
matrix_k <- as.matrix(t(heat_map_k))  # transpose
heatmap(matrix_k, Rowv = NA, Colv = NA, scale = "none", main = "k1 vs. k2", xlab = "k1", ylab = "k2", col = jOrRdPalette)
mtext("Constants", 1, line = -12, adj = 1)
mtext("vac1 = vac2 = .5", 1, line = -11, adj = 1)
mtext("s0 = v1 = v2 = 3300", 1, line = -10, adj = 1)

############################
# Compare vaccine efficacy #
############################
par <- par_original # reset parameters
par$k1 = 0.5  # set the relative infectiousness of once-vax and twice-vax to the same
par$k2 = 0.5

heat_map_vac <- data.frame(matrix(rep(0, bins*bins), nrow=bins))

for (x in 1:bins){             # changes in vac1
  par$vac2 = .1*x
  rownames(heat_map_vac)[x] <- as.character(par$vac2)
  for (y in 1:x){              # changes in vac2
    par$vac1 = .1*y
    colnames(heat_map_vac)[y] <- as.character(par$vac1)
    for (z in 1:repeats){      # run trials and average number recovered
      trial <- Run_Model_S(par=par, inits = c(3300,3300,3300,0,0,0,40,30,30,0,0,0), campaign=0)
      heat_map_vac[x,y] <- heat_map_vac[x,y] + sum(trial[nrow(trial), c("R_0", "R_1","R_2")])
    }
  }
}
heat_map_vac <- heat_map_vac/repeats

# create a heatmap 
matrix_vac <- as.matrix(t(heat_map_vac))
heatmap(matrix_vac, Rowv = NA, Colv = NA, scale = "none", main = "vac2 vs. vac1", xlab = "vac2", ylab = "vac1", col = jOrRdPalette)
mtext("Constants", 1, line = -12, adj = 1)
mtext("k1 = k2 = .5", 1, line = -11, adj = 1)
mtext("s0 = v1 = v2 = 3300", 1, line = -10, adj = 1)

###################################################
# Compare relative number preemptively vaccinated #
###################################################
par <- par_original # reset parameters

heat_map_v <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
colnames(heat_map_v) <- rep(NA,ncol(heat_map_v))
rownames(heat_map_v) <- as.character((1:10)/100) #just a work-around to initialize

for (x in 1:bins){             # changes in number of vaccines given to first dose
  v1 = x*100
  rownames(heat_map_v)[x] <- as.character(v1)
  for (y in 1:bins){              # changes in number of vaccines given in second dose
    v2 = y*100
    s0 = par$ntot  - 100 - v1 - v2
    colnames(heat_map_v)[y] <- as.character(v2)
    for (z in 1:repeats){      # run trials and average number recovered
      trial <- Run_Model_S(par=par, inits = c(s0,v1,v2,0,0,0,40,30,30,0,0,0), campaign=0)
      heat_map_v[x,y] <- heat_map_v[x,y] + sum(trial[nrow(trial), c("R_0", "R_1","R_2")])
    }
  }
}

heat_map_v <- heat_map_v/repeats

# create a heatmap 
matrix_v <- as.matrix(t(heat_map_v))
heatmap(matrix_v, Rowv = NA, Colv = NA, scale = "none", main = "v1 vs. v2", xlab = "v1", ylab = "v2", col = jOrRdPalette)
mtext("Constants", 1, line = -12, adj = 1)
mtext("k1 = .2, k2 = .5", 1, line = -11, adj = 1)
mtext("vac1 = .5, vac2 = .8", 1, line = -10, adj = 1)

##############################################################
# Compare relative infectiousness and preemptive vaccination #
##############################################################
par <- par_original # reset parameters

heat_map_k_v <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
rownames(heat_map_k_v) <- as.character((1:10)/100) #just a work-around to initialize
par$vac1 = 0.5  # set the protective efficacy of either dose to 0.5
par$vac = 0.5
par$k2 = .9
for (x in 1:bins){
  par$k1 = .1*x
  ratio_k = par$k1/par$k2
  rownames(heat_map_k_v)[x] <- paste("[",par$k1,",",par$k2,"]")
  for (y in 1:bins){
    v1 = 5000 - 500*(y-1)
    v2 = (5000 - v1)/2
    ratio_v = v1/v2
    s0 = par$ntot - 100 - v1 - v2
    colnames(heat_map_k_v)[y] <- paste("[",v1,",",v2,"]")
    for (z in 1:repeats){      # run trials and average number recovered
      trial <- Run_Model_S(par=par, inits = c(s0,v1,v2,0,0,0,40,30,30,0,0,0), campaign=0)
      heat_map_k_v[x,y] <- heat_map_k_v[x,y] + sum(trial[nrow(trial), c("R_0", "R_1","R_2")])
    }
  }
}
heat_map_k_v <- heat_map_k_v/repeats

# create a heatmap 
matrix_k_v <- as.matrix(t(heat_map_k_v))
heatmap(matrix_k_v, Rowv = NA, Colv = NA, scale = "none", main = "k1/k2 vs. v1/v2", col = jOrRdPalette)
mtext("Constants", 1, line = -12, adj = 1.05)
mtext("vac1 = vac2 = .5", 1, line = -11, adj = 1.05)
mtext("[k1,k2]", 2, line = -4.5, adj = -.7)
mtext("[v1,v2]", 1, line = -12.7, adj = .72)
#mtext("k2 = .1", 1, line = -10, adj = 1)
# Bottom left corner is if you give all your doses to v1 and if you are just as infectious if you got vaccinated once as twice
# Top right corner is if you twice-vaccinate 4.5 times more people than once-vaccinate and if your reduction in infectiousness is 10 times greater after the second dose

################################################################
# Compare relative vaccine efficacy and preemptive vaccination #
###############################################################
par <- par_original # reset parameters

heat_map_vac_v <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
rownames(heat_map_vac_v) <- as.character((1:10)/100) #just a work-around to initialize
par$k1 = 0.5    # set relative infectiousness to 0.5 for each
par$k2 = 0.5
par$vac1 = .1
for (x in 1:bins){
  par$vac2 = .1*x
  ratio_vac = par$vac1/par$vac2
  rownames(heat_map_vac_v)[x] <- paste("[",par$vac1,",",par$vac2,"]")  # calculating the row names directly
  for (y in 1:bins){
    v1 = 5000 - 500*(y-1)
    v2 = (5000 - v1)/2
    ratio_v = v1/v2
    colnames(heat_map_vac_v)[y] <- paste("[",v1,",",v2,"]")   # calculating the column names directly
    s0 = par$ntot - 100 - v1 - v2
    for (z in 1:repeats){      # run trials and average number recovered
      trial <- Run_Model_S(par=par, inits = c(s0,v1,v2,0,0,0,40,30,30,0,0,0), campaign=0)
      heat_map_vac_v[x,y] <- heat_map_vac_v[x,y] + sum(trial[nrow(trial), c("R_0", "R_1","R_2")])
    }
  }
}
heat_map_vac_v <- heat_map_vac_v/repeats

# create a heatmap 
matrix_vac_v <- as.matrix(t(heat_map_vac_v))
heatmap(matrix_vac_v, Rowv = NA, Colv = NA, scale = "none", main = "vac1/vac2 vs. v1/v2", col = jOrRdPalette)
mtext("Constants", 1, line = -12, adj = 1)
mtext("k1 = k2 = .5", 1, line = -11, adj = 1)
#mtext("vac1 = .1", 1, line = -10, adj = 1)
mtext("[vac1,vac2]", 2, line = -4.5, adj = -1.05)
mtext("[v1,v2]", 1, line = -12.7, adj = .72)
# Bottom left corner is if we give all our doses to the first dose and the doses are equally effectve
# Top right corner is if we create 4.5 times more twice-vax people than once-vax people and if the first dose is 10% effective and the second dose is 100% effective
# Note the bottom row should really be equal since the efficacy of the first dose is constant The top row and right column should have the fastest change

########################################################
# Compare relative vaccine efficacy and infectiousness #
########################################################
par <- par_original # reset parameters

heat_map_vac_k <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
rownames(heat_map_vac_k) <- as.character((1:10)/100) #just a work-around to initialize
par$vac1 = .1
for (x in 1:bins){
  par$vac2 = .1*x
  ratio_vac = par$vac1/par$vac2
  rownames(heat_map_vac_k)[x] <- paste("[",par$vac1,",",par$vac2,"]")
  par$k2 = .9
  for (y in 1:bins){
    par$k1 = .1*y
    ratio_k = par$k1/par$k2
    colnames(heat_map_vac_k)[y] <- paste("[",par$k1,",",par$k2,"]")
    for (z in 1:repeats){      # run trials and average number recovered
      trial <- Run_Model_S(par=par, inits = c(3300,3300,3300,0,0,0,40,30,30,0,0,0), campaign=0)
      heat_map_vac_k[x,y] <- heat_map_vac_k[x,y] + sum(trial[nrow(trial), c("R_0", "R_1","R_2")])
    }
  }
}
heat_map_vac_k <- heat_map_vac_k/repeats

# create a heatmap 
matrix_vac_k <- as.matrix(t(heat_map_vac_k))
heatmap(matrix_vac_k, Rowv = NA, Colv = NA, scale = "none", main = "vac1/vac2 vs. k1/k2", col = jOrRdPalette)
mtext("Constants", 1, line = -12, adj = 1)
#mtext("vac1 = .1", 1, line = -11, adj = 1)
#mtext("k2 = .9", 1, line = -10, adj = 1)
mtext("[vac1,vac2]", 2, line = -4.5, adj = -1.05)
mtext("[k1,k2]", 1, line = -12.7, adj = .72)
mtext("s0 = v1 = v2 = 3300", 1, line = -11, adj = 1)

# The bottom left corner  is if the first dose is 10% effective, the second dose is 100% effective, and both reduce infectiousness by 90%
# The bottom right corner is if both doses are 10% effective and both reduce infectiousness by 90%
# The top right corner    is if both doses are 10% effective, the first dose does not reduce infectiousness but the second dose reduces it by 90%
# The top left corner     is if the first dose is 10% effective, the second dose is 100% effective, the first dose does not reduce infectiousness but the second dose reduces it by 90%
# Very interesting split in the data around k1/k2=4. 

## Note: Reactive vaccinations are no longer useful comparisons because
## we will choose a 1 or 2 dose campaign, not a mixture; earlier code has
## been modified to reflect this; changes must be made for reactive
## comparison code to run

########################
# Reactive Vaccination #
########################
par <- par_original # reset parameters

par$v1_day = 10
heat_map_reac <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
rownames(heat_map_reac) <- as.character((1:10)/100) #just a work-around to initialize
for (x in 1:bins){
  par$v1_count = 100*x
  rownames(heat_map_reac)[x] <- as.character(par$v1_count)
  for (y in 1:bins){
    par$v2_count = 100*y
    colnames(heat_map_reac)[y] <- as.character(par$v2_count)
    for (z in 1:repeats){      # run trials and average number recovered
      trial <- Run_Model_S(par=par, inits = c(9990,0,0,0,0,0,10,0,0,0,0,0), campaign=2)
      heat_map_reac[x,y] <- heat_map_reac[x,y] + sum(trial[nrow(trial), c("R_0", "R_1","R_2")])
    }
  }
}
heat_map_reac <- heat_map_reac/repeats

# create a heatmap 
matrix_reac <- as.matrix(t(heat_map_reac))
heatmap(matrix_reac, Rowv = NA, Colv = NA, scale = "none", main = "v1_count vs. v2_count", xlab = "v1_count", ylab = "v2_count", col = jOrRdPalette)
mtext("Constants", 1, line = -12, adj = 1)
mtext("vac1 = .5, vac2 = .8", 1, line = -11, adj = 1)
mtext("k1 = .2, k2 = .5", 1, line = -10, adj = 1)
mtext("s0 = 9990", 1, line = -9, adj = 1)

############################################################
# Compare relative infectiousness and reactive vaccination #
############################################################
par <- par_original # reset parameters
par$v1_day = 10

heat_map_k_reac <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
rownames(heat_map_k_reac) <- as.character((1:10)/100) #just a work-around to initialize
par$vac1 = 0.5  # set the protective efficacy of either dose to 0.5
par$vac = 0.5
par$k2 = .9
for (x in 1:bins){
  par$k1 = .1*x
  ratio_k = par$k1/par$k2
  rownames(heat_map_k_reac)[x] <- paste("[",par$k1,",",par$k2,"]")
  for (y in 1:bins){
    par$v1_count = 5000 - 500*(y-1)
    par$v2_count = (5000 - par$v1_count)/2
    colnames(heat_map_k_reac)[y] <- paste("[",par$v1_count,",",par$v2_count,"]")
    for (z in 1:repeats){      # run trials and average number recovered
      trial <- Run_Model_S(par=par, inits = c(9999,0,0,0,0,0,1,0,0,0,0,0), campaign=2)
      heat_map_k_reac[x,y] <- heat_map_k_reac[x,y] + sum(trial[nrow(trial), c("R_0", "R_1","R_2")])
    }
  }
}
heat_map_k_reac <- heat_map_k_reac/repeats

# create a heatmap 
matrix_k_reac <- as.matrix(t(heat_map_k_reac))
heatmap(matrix_k_reac, Rowv = NA, Colv = NA, scale = "none", main = "k vs. v_count", col = jOrRdPalette)
mtext("Constants", 1, line = -12, adj = 1.05)
mtext("vac1 = vac2 = .5", 1, line = -11, adj = 1.05)
mtext("s0 = 9900", 1, line = -10, adj = 1)
mtext("[k1,k2]", 2, line = -4.5, adj = -.7)
mtext("[v1_count,v2_count]", 1, line = -12.7, adj = .82)

##############################################################
# Compare relative vaccine efficacy and reactive vaccination #
##############################################################
par <- par_original # reset parameters
par$v1_day = 5

heat_map_vac_reac <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
rownames(heat_map_vac_reac) <- as.character((1:10)/100) #just a work-around to initialize
par$k1 = 0.5    # set relative infectiousness to 0.5 for each
par$k2 = 0.5
par$vac1 = .1
for (x in 1:bins){
  par$vac2 = .1*x
  rownames(heat_map_vac_reac)[x] <- paste("[",par$vac1,",",par$vac2,"]")  # calculating the row names directly
  for (y in 1:bins){
    par$v1_count = 5000 - 500*(y-1)
    par$v2_count = (5000 - par$v1_count)/2
    colnames(heat_map_vac_reac)[y] <- paste("[",par$v1_count,",",par$v2_count,"]")   # calculating the column names directly
    for (z in 1:repeats){      # run trials and average number recovered
      trial <- Run_Model_S(par=par, inits = c(9999,0,0,0,0,0,1,0,0,0,0,0), campaign=2)
      heat_map_vac_reac[x,y] <- heat_map_vac_reac[x,y] + sum(trial[nrow(trial), c("R_0", "R_1","R_2")])
    }
  }
}
heat_map_vac_reac <- heat_map_vac_reac/repeats

# create a heatmap 
matrix_vac_reac <- as.matrix(t(heat_map_vac_reac))
heatmap(matrix_vac_reac, Rowv = NA, Colv = NA, scale = "none", main = "vac vs. v_count", col = jOrRdPalette)
mtext("Constants", 1, line = -12, adj = 1)
mtext("k1 = k2 = .5", 1, line = -11, adj = 1)
mtext("s0 = 9900", 1, line = -10, adj = 1)
mtext("[vac1,vac2]", 2, line = -4.5, adj = -1.1)
mtext("[v1_count,v2_count]", 1, line = -12.7, adj = .82)

###########################################
# Compare delay time and vaccine efficacy #
###########################################
par <- par_original # reset parameters
par$v1_day = 10
par$vac1 = .1

heat_map_1dose <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
heat_map_2dose <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
heat_map_0dose <- data.frame(c(rep(0, bins)))
heat_map_2dose_half <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
rownames(heat_map_1dose) <- as.character((1:10)*10)
rownames(heat_map_2dose) <- as.character((1:10)*10)
rownames(heat_map_0dose) <- as.character((1:10)*100)
rownames(heat_map_2dose_half) <- as.character((1:10)*10)

for (z in 1:repeats){      # run trials and average number recovered
  for (x in 1:bins){
    par$v_delay = 10*x
    rownames(heat_map_1dose)[x] <- par$v_delay
    rownames(heat_map_2dose)[x] <- par$v_delay
    rownames(heat_map_0dose)[x] <- par$v_delay
    rownames(heat_map_2dose_half)[x] <- par$v_delay
    for (y in 1:bins){
      par$vac2 = .1*y
      colnames(heat_map_1dose)[y] <- paste("[",par$vac1,",",par$vac2,"]")
      colnames(heat_map_2dose)[y] <- paste("[",par$vac1,",",par$vac2,"]")
      colnames(heat_map_2dose_half)[y] <- paste("[",par$vac1,",",par$vac2,"]")
      
      par$v_count = 2000
      trial1 <- Run_Model_S(par=par, inits = c(9990,0,0,0,0,0,10,0,0,0,0,0), campaign=1)
      heat_map_1dose[x,y] <- heat_map_1dose[x,y] + sum(trial1[nrow(trial1), c("R_0", "R_1","R_2")])
      trial2 <- Run_Model_S(par=par, inits = c(9990,0,0,10,0,0,0,0,0), campaign=2)
      heat_map_2dose[x,y] <- heat_map_2dose[x,y] + sum(trial2[nrow(trial2), c("R_0", "R_1","R_2")])
      
      par$v_count = par$v_count/2
      trial2_half <- Run_Model_S(par=par, inits = c(9990,0,0,0,0,0,10,0,0,0,0,0), campaign=2)
      heat_map_2dose_half[x,y] <- heat_map_2dose_half[x,y] + sum(trial2_half[nrow(trial2_half), c("R_0", "R_1","R_2")])
    }
    par$v_count = par$v_count*2
    trial0 <- Run_Model_S(par=par, inits = c(9990,0,0,0,0,0,10,0,0,0,0,0), campaign=0)
    heat_map_0dose[x,1] <- heat_map_0dose[x,1] + sum(trial0[nrow(trial0), c("R_0", "R_1", "R_2")])
  }
}

heat_map_1dose <- heat_map_1dose/repeats
heat_map_2dose <- heat_map_2dose/repeats
heat_map_0dose <- heat_map_0dose/repeats
heat_map_2dose_half <- heat_map_2dose_half/repeats

# create comparisons
matrix_1dose_both <- 2*as.matrix(t(heat_map_1dose))
matrix_2dose <- as.matrix(t(heat_map_2dose))
mat0 <- as.matrix(t(heat_map_0dose))
matrix_0dose <- rbind(mat0,mat0,mat0,mat0,mat0,mat0,mat0,mat0,mat0,mat0)
matrix_2dose_one <- matrix_2dose + matrix_0dose
matrix_2dosehalf_both <- 2*as.matrix(t(heat_map_2dose_half))
matrix_comparison <- matrix(, nrow = bins, ncol = 10)
for (x in 1:bins){
  for (y in 1:bins){
    if (min(matrix_1dose_both[x,y],matrix_2dose_one[x,y],matrix_2dosehalf_both[x,y]) == matrix_1dose_both[x,y]){
      matrix_comparison[x,y] = 1
    }
    else if (min(matrix_1dose_both[x,y],matrix_2dose_one[x,y],matrix_2dosehalf_both[x,y]) == matrix_2dose_one[x,y]){
      matrix_comparison[x,y] = 2
    }
    else if (min(matrix_1dose_both[x,y],matrix_2dose_one[x,y],matrix_2dosehalf_both[x,y]) == matrix_2dosehalf_both[x,y]){
      matrix_comparison[x,y] = 3
    }
    else{
      matrix_comparison[x,y] = NA
    }
  }
}
rownames(matrix_comparison) <- colnames(heat_map_1dose)
colnames(matrix_comparison) <- rownames(heat_map_1dose)
# create a heatmap 
heatmap(matrix_comparison, Rowv = NA, Colv = NA, scale = "none", main = "v_delay vs. vaccine efficacy", col = jOrRdPalette)
# mtext("Constants", 1, line = -12, adj = 1.05)
# mtext("vac1 = .5, vac2 = .8", 1, line = -11, adj = 1.05)
# mtext("s0 = 9900", 1, line = -10, adj = 1)
mtext("[v_delay]", 2, line = -4.5, adj = -.7)
mtext("[vac1,vac2]", 1, line = -12.7, adj = .82)

#############################
# Probability of Extinction #
#############################
run_time <- proc.time()
Pr_extinct <- matrix(rep(0,300),nrow=100)
#assume we can vaccinate almost everyone
inits <- c(9999,0,0,0,0,0,1,0,0,0,0,0)
par=par_original
doses = 8000
for (h in 1:10){
  cutoff <- 1000*h
  for (i in 1:50){
    for (j in 1:10){
      par$v_count = doses
      if (j == 1){
        #preemptively vaccinate
        par$v1_day = 1
        par$v_delay = 1
        trial1 <- Run_Model_S(par=par, inits=inits, campaign=1)
        trial2 <- Run_Model_S(par=par, inits=inits, campaign=2)
        par$v_count = 2*par$v_count
        trial2_double <- Run_Model_S(par=par, inits=inits, campaign=2)
      }
      else{
        #normal
        par$v1_day = 5*(j-1)
        par$v_delay = 10
        trial1 <- Run_Model_S(par=par, inits=inits, campaign=1)
        trial2 <- Run_Model_S(par=par, inits=inits, campaign=2)
        par$v_count = 2*par$v_count
        trial2_double <- Run_Model_S(par=par, inits=inits, campaign=2)
      }
      if (sum(trial1[nrow(trial1), c("R_0", "R_1", "R_2")]) < cutoff){
        Pr_extinct[(h-1)*10+j,1] <- Pr_extinct[(h-1)*10+j,1] + 1
      }
      if (sum(trial2[nrow(trial2), c("R_0", "R_1", "R_2")]) < cutoff){
        Pr_extinct[(h-1)*10+j,2] <- Pr_extinct[(h-1)*10+j,2] + 1
      }
      if (sum(trial2_double[nrow(trial2_double), c("R_0", "R_1", "R_2")]) < cutoff){
        Pr_extinct[(h-1)*10+j,3] <- Pr_extinct[(h-1)*10+j,3] + 1
      }
    }
  }
}
Pr_extinct <- Pr_extinct/i
View(Pr_extinct)
colnames(Pr_extinct) <- c("1 dose", "2 dose", "2 dose double")

run_time <- proc.time() - run_time
run_time

# determine marginal probability of extinction of 1 dose or 2 dose double campaign
Pr_extinct_marg <- Pr_extinct
Pr_extinct_marg[,1] <- Pr_extinct[,2]
Pr_extinct_marg[,2] <- Pr_extinct[,1] - Pr_extinct[,2]
Pr_extinct_marg[,3] <- Pr_extinct[,3] - Pr_extinct_marg[,1] - Pr_extinct_marg[,2]
colnames(Pr_extinct_marg) <- c("2 dose", "1 dose", "2 dose double")
rownames(Pr_extinct_marg) <- c(rep(5*(1:10) - 5,10))

# make stacked bar graphs
stacked_bar1_gg <- melt(Pr_extinct_marg[1:10,])
ggplot(stacked_bar1_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(1000) + labs(x="v1_day", y="Probability of Extinction")
stacked_bar2_gg <- melt(Pr_extinct_marg[11:20,])
ggplot(stacked_bar2_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(2000) + labs(x="v1_day", y="Probability of Extinction")
stacked_bar3_gg <- melt(Pr_extinct_marg[21:30,])
ggplot(stacked_bar3_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(3000) + labs(x="v1_day", y="Probability of Extinction")
stacked_bar4_gg <- melt(Pr_extinct_marg[31:40,])
ggplot(stacked_bar4_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(4000) + labs(x="v1_day", y="Probability of Extinction")
stacked_bar5_gg <- melt(Pr_extinct_marg[41:50,])
ggplot(stacked_bar5_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(5000) + labs(x="v1_day", y="Probability of Extinction")
stacked_bar6_gg <- melt(Pr_extinct_marg[51:60,])
ggplot(stacked_bar6_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(6000) + labs(x="v1_day", y="Probability of Extinction")
stacked_bar7_gg <- melt(Pr_extinct_marg[61:70,])
ggplot(stacked_bar7_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(7000) + labs(x="v1_day", y="Probability of Extinction")
stacked_bar8_gg <- melt(Pr_extinct_marg[71:80,])
ggplot(stacked_bar8_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(8000) + labs(x="v1_day", y="Probability of Extinction")
stacked_bar9_gg <- melt(Pr_extinct_marg[81:90,])
ggplot(stacked_bar9_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(9000) + labs(x="v1_day", y="Probability of Extinction")
stacked_bar10_gg <- melt(Pr_extinct_marg[91:100,])
ggplot(stacked_bar10_gg, aes(x=factor(Var1), y=value, fill=factor(Var2))) + geom_bar(stat="identity") + ggtitle(10000) + labs(x="v1_day", y="Probability of Extinction")

#####################
# Epidemic Duration #
#####################
inits_s <- c(24999,0,0,0,0,0,1,0,0,0,0,0)
trial_size <- 100
vacs <- 1
duration_mat <- data.frame(c(rep(0,trial_size*(vacs+1))))
par <- par_original
par$v1_day <- 10
par$v_delay <- 1
par1 <- par
par1$v_delay <- 10
par2 <- par
par2$v_delay <- 20
for (i in 1:trial_size){
  trial <- Run_Model_S(par = par, inits = inits_s, campaign = 1)
  j <- 1
  while(j<nrow(trial)&&(trial[j,5]!=0||trial[j,6]!=0||trial[j,7]!=0||trial[j,8]!=0||trial[j,9]!=0||trial[j,10]!=0)){
    j = j + 1
  }
  duration_mat[i,1] <- trial[j,1]
  duration_mat[i,2] <- "1 day"
  
  trial_1 <- Run_Model_S(par = par1, inits = inits_s, campaign = 1)
  k <- 1
  while(k<nrow(trial_1)&&(trial_1[k,5]!=0||trial_1[k,6]!=0||trial_1[k,7]!=0||trial_1[k,8]!=0||trial_1[k,9]!=0||trial_1[k,10]!=0)){
    k = k + 1
  }
  duration_mat[i+trial_size,1] <- trial_1[k,1]
  duration_mat[i+trial_size,2] <- "10 day"
  
  trial_2 <- Run_Model_S(par = par2, inits = inits_s, campaign = 2)
  l <- 1
  while(l<nrow(trial_2)&&(trial_2[l,5]!=0||trial_2[l,6]!=0||trial_2[l,7]!=0||trial_2[l,8]!=0||trial_2[l,9]!=0||trial_2[l,10]!=0)){
    l = l + 1
  }
  duration_mat[i+2*trial_size,1] <- trial_2[l,1]
  duration_mat[i+2*trial_size,2] <- "20 day"
}
colnames(duration_mat) <- c("duration", "vaccines")
#qplot(duration,data=duration_mat, geom="histogram", binwidth=1, fill="blue")
h <- ggplot(duration_mat, aes(x=duration))
h + geom_histogram(binwidth=2, fill="red", colour="black") + facet_grid(vaccines ~ .) 

############################
# Cumulative Epidemic Size #
############################
parms <- parms_original
t_final <- 200
inits <- c(S_0=24999,
           S_1=0,
           S_2=0,
           E_0=0,
           E_1=0,
           E_2=0,
           I_0=1,
           I_1=0,
           I_2=0,
           R_0=0,
           R_1=0,
           R_2=0,
           N  =25000)
start <- 6
delay <- 14 #Note, Corey changed this to 14 on Dec 14
end <- t_final - delay - 1
time_steps <- end - start + 1
epidemic_size <- data.frame(matrix(rep(0,time_steps*4),nrow=time_steps))
epidemic_size[,1] <- start:end
colnames(epidemic_size) <- c("time", "AB_1dose", "A_2dose", "AB_2dosehalf")
# assume a fixed number of vaccines = enough to vaccinate everyone once
for (i in start:end){
  parms["v1_day"] <- i
  parms["v2_day"] <- parms["v1_day"] + delay
  # 1 dose AB
  parms["v1_count"] <- inits["N"]
  parms["v2_count"] <- 0
  simulation0 <- Run_Model_D(inits, dt, parms=parms)
  simulation0.df <- as.data.frame(simulation0)
  simulation0.df$D <- 1-apply(simulation0.df[,2:length(simulation0.df)], 1, sum)
  epidemic_size[i,2] <- 2*sum(simulation0.df[nrow(simulation0.df), c("R_0", "R_1","R_2")])
  # 2 dose A
  parms["v1_count"] <- inits["N"]
  parms["v2_count"] <- inits["N"]
  simulation1a <- Run_Model_D(inits, dt, parms=parms)
  simulation1a.df <- as.data.frame(simulation1a)
  simulation1a.df$D <- 1-apply(simulation1a.df[,2:length(simulation1a.df)], 1, sum)
  epidemic_size[i,3] <- sum(simulation1a.df[nrow(simulation1a.df), c("R_0", "R_1","R_2")]) 
  # 2 dose half AB
  parms["v1_count"] <- inits["N"]/2
  parms["v2_count"] <- inits["N"]/2
  simulation2 <- Run_Model_D(inits, dt, parms=parms)
  simulation2.df <- as.data.frame(simulation2)
  simulation2.df$D <- 1-apply(simulation2.df[,2:length(simulation2.df)], 1, sum)
  epidemic_size[i,4] <- 2*sum(simulation2.df[nrow(simulation2.df), c("R_0", "R_1","R_2")])
}
# 0 dose (same for all v1_days and v2_days)
parms["v1_count"] <- 0
parms["v2_count"] <- 0
simulation1b <- Run_Model_D(inits, dt, parms=parms)
simulation1b.df <- as.data.frame(simulation1b)
simulation1b.df$D <- 1-apply(simulation1b.df[,2:length(simulation1b.df)], 1, sum)
epidemic_size[,3] <- epidemic_size[,3] + sum(simulation1b.df[nrow(simulation1b.df), c("R_0", "R_1","R_2")])
ggplot(data=epidemic_size, aes(x=time)) + 
  geom_line(aes(x=time, y=AB_1dose), col="blue") + 
  geom_line(aes(x=time, y=A_2dose), col="red") + 
  geom_line(aes(x=time, y=AB_2dosehalf), col="green") + 
  theme_bw() + xlab("Time Until First Vaccine") + ylab("Cumulative Epidemic Size") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Epidemic Size According to Campaign")

#Corey, Dec 14#
parms <- parms_original
parms["v1_count"] <- 0
parms["v2_count"] <- 0
simulation <- Run_Model_D(inits, dt, parms=parms)
simulation.df <- as.data.frame(simulation)
simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
Plot_SEIR(simulation.df)
ggplot(data=epidemic_size, aes(x=time)) + 
  geom_bar(data=simulation.df[1:nrow(epidemic_size),], aes(x=time, y=I_0*2), stat="identity", fill="lightgrey", color="lightgrey") +
  geom_line(aes(x=time, y=AB_1dose), col="blue") + 
  geom_line(aes(x=time, y=A_2dose), col="red") + 
  geom_line(aes(x=time, y=AB_2dosehalf), col="green") + 
  theme_bw() + xlab("Time Until First Vaccine") + ylab("Cumulative Epidemic Size") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Epidemic Size According to Campaign")

########################################################
# Choice Heat Map: Vaccine Timing and Vaccine Efficacy #
########################################################
bins <- 11
vaccines1 <- 3000
parms <- parms_original
# different inits for different preemptive vaccinations
inits_a <- c(S_0=max(9999-vaccines1/2,0),
             S_1=min(9999,vaccines1/2),
             S_2=0,
             E_0=0,
             E_1=0,
             E_2=0,
             I_0=1,
             I_1=0,
             I_2=0,
             R_0=0,
             R_1=0,
             R_2=0,
             N  =10000)
inits_b <- c(S_0=max(9999-vaccines1/2,0),
           S_1=0,
           S_2=min(9999,vaccines1/2),
           E_0=0,
           E_1=0,
           E_2=0,
           I_0=1,
           I_1=0,
           I_2=0,
           R_0=0,
           R_1=0,
           R_2=0,
           N  =10000)
inits_c <- c(S_0=max(9999-vaccines1/4,0),
             S_1=0,
             S_2=min(9999,vaccines1/4),
             E_0=0,
             E_1=0,
             E_2=0,
             I_0=1,
             I_1=0,
             I_2=0,
             R_0=0,
             R_1=0,
             R_2=0,
             N  =10000)
inits_reactive <- c(S_0=9999,
           S_1=0,
           S_2=0,
           E_0=0,
           E_1=0,
           E_2=0,
           I_0=1,
           I_1=0,
           I_2=0,
           R_0=0,
           R_1=0,
           R_2=0,
           N  =10000)

heat_map_AB_1dose <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
heat_map_A_2dose <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
heat_map_AB_2dosehalf <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
rownames(heat_map_AB_1dose) <- c(1:bins)
colnames(heat_map_AB_1dose) <- c(1:bins)
for (x in 1:bins){
  rownames(heat_map_AB_1dose)[x] <- x/10
  for (y in 1:bins){
    colnames(heat_map_AB_1dose)[y] <- y/100
  }
}
for (x in 1:bins){
  parms["vac1"] = (parms["vac2"]*(x-1))/(bins-1)
  rownames(heat_map_AB_1dose)[x] <- as.character(parms["vac1"])
  for (y in 1:bins){
    if (y == 1){
      # preemptively vaccinate
      inits_0 <- inits_a
      inits_1 <- inits_b
      inits_2 <- inits_c
      inits_00 <- inits_reactive
      vaccines <- 0 # no more available vaccines
      parms["v1_day"] <- 6 #these numbers are placeholders
      parms["v2_day"] <- 8
      colnames(heat_map_AB_1dose)[y] <- "preemptive"
    } else{
      vaccines <- vaccines1 # reset number of vaccines
      inits_0 <- inits_reactive
      inits_1 <- inits_reactive
      inits_2 <- inits_reactive
      inits_00 <- inits_reactive
      parms["v1_day"] <- round(y/bins*88,0)
      parms["v2_day"] <- parms["v1_day"] + 14
      colnames(heat_map_AB_1dose)[y] <- as.character(parms["v1_day"])
    }
    # 1 dose AB
    parms["v1_count"] <- vaccines/2
    parms["v2_count"] <- 0
    simulation0 <- Run_Model_D(inits_0, dt, parms=parms)
    simulation0.df <- as.data.frame(simulation0)
    simulation0.df$D <- 1-apply(simulation0.df[,2:length(simulation0.df)], 1, sum)
    heat_map_AB_1dose[x,y] <- sum(simulation0.df[nrow(simulation0.df), c("R_0", "R_1","R_2")])
    # 2 dose A
    parms["v1_count"] <- vaccines/2
    parms["v2_count"] <- vaccines/2
    simulation1a <- Run_Model_D(inits_1, dt, parms=parms)
    simulation1a.df <- as.data.frame(simulation1a)
    simulation1a.df$D <- 1-apply(simulation1a.df[,2:length(simulation1a.df)], 1, sum)
    heat_map_A_2dose[x,y] <- sum(simulation1a.df[nrow(simulation1a.df), c("R_0", "R_1","R_2")])
    # 2 dose half AB
    parms["v1_count"] <- vaccines/4
    parms["v2_count"] <- vaccines/4
    simulation2 <- Run_Model_D(inits_2, dt, parms=parms)
    simulation2.df <- as.data.frame(simulation2)
    simulation2.df$D <- 1-apply(simulation2.df[,2:length(simulation2.df)], 1, sum)
    heat_map_AB_2dosehalf[x,y] <- sum(simulation2.df[nrow(simulation2.df), c("R_0", "R_1","R_2")])
  }
}
# 0 dose B
parms["v1_count"] <- 0
parms["v2_count"] <- 0
simulation1b <- Run_Model_D(inits_0, dt, parms=parms)
simulation1b.df <- as.data.frame(simulation1b)
simulation1b.df$D <- 1-apply(simulation1b.df[,2:length(simulation1b.df)], 1, sum)
heat_map_B_0dose <- matrix(rep(sum(simulation1b.df[nrow(simulation1b.df), c("R_0", "R_1","R_2")]),bins*bins), nrow=bins)



# create comparisons
matrix_1 <- 2*as.matrix((heat_map_AB_1dose))
matrix_2_A <- as.matrix((heat_map_A_2dose))
matrix_0_B <- as.matrix((heat_map_B_0dose))
matrix_2 <- matrix_2_A + matrix_0_B
matrix_2.5 <- 2*as.matrix((heat_map_AB_2dosehalf))
matrix_choice <- matrix(, nrow = bins, ncol = bins)
colnames(matrix_choice) <- colnames(heat_map_AB_1dose)
rownames(matrix_choice) <- rownames(heat_map_AB_1dose)
for (v in 1:bins){
  for (w in 1:bins){
    if (min(matrix_1[v,w],matrix_2[v,w],matrix_2.5[v,w]) == matrix_1[v,w]){
      matrix_choice[v,w] = 1
    }
    else if (min(matrix_1[v,w],matrix_2[v,w],matrix_2.5[v,w]) == matrix_2[v,w]){
      matrix_choice[v,w] = 2
    }
    else if (min(matrix_1[v,w],matrix_2[v,w],matrix_2.5[v,w]) == matrix_2.5[v,w]){
      matrix_choice[v,w] = 3
    }
    else{
      matrix_choice[v,w] = 4
    }
  }
}
# create a heatmap 
myCol <- c("blue", "red", "green", "black")
## Defining breaks for the color scale
myBreaks <- c(.5, 1.5, 2.5, 3.5, 4.5)

graph_title <- "Vaccine Timing vs. Vaccine Efficacy"
heatmap(matrix_choice, scale="none", Rowv=NA, Colv=NA, main = graph_title, col = myCol, breaks = myBreaks, xlab="v1_day", ylab="vac1") 

mtext("Constants", 1, line = -13, adj = 1)
mtext("k1 = k2 = 0", 1, line = -12, adj = 1)
mtext("vac2 = .8", 1, line = -11, adj = 1)
text <- paste(vaccines1, " vaccines")
mtext(text, 1, line = -10, adj = 1)

matrix_1_vs_2.5 <- matrix_2.5 - matrix_1
colnames(matrix_1_vs_2.5) <- colnames(heat_map_AB_1dose)
rownames(matrix_1_vs_2.5) <- rownames(heat_map_AB_1dose)
#heatmap(matrix_1_vs_2.5, scale="none", Rowv=NA, Colv=NA, main="Additional People Saved, 1 vs. 2.5", col = jPRGnPalette, xlab="v1_day", ylab="vac1")
#View(matrix_1_vs_2.5)

# matrix_1_vs_2 <- matrix_2 - matrix_1
# colnames(matrix_1_vs_2) <- colnames(heat_map_AB_1dose)
# rownames(matrix_1_vs_2) <- rownames(heat_map_AB_1dose)
# heatmap(matrix_1_vs_2, scale="none", Rowv=NA, Colv=NA, main="Additional People Saved, 1 vs. 2.5", col = jRdBuPalette, xlab="k1", ylab="vac1")

scaling_factor <- abs(max(matrix_1_vs_2.5)/min(matrix_1_vs_2.5))
matrix_scaled <- matrix_1_vs_2.5
colnames(matrix_scaled) <- colnames(heat_map_AB_1dose)
rownames(matrix_scaled) <- rownames(heat_map_AB_1dose)
graph_title2 <- "Additional People Saved, 1 vs. 2.5"
max_saved_green <- paste("max saved green = ", round(-min(matrix_1_vs_2.5)))
vacs <- paste(vaccines, " vaccines")
if (max(matrix_1_vs_2.5) > 0){
  for (x in 1:bins){
    for (y in 1:bins){
      if (matrix_1_vs_2.5[x,y] < 0){
        matrix_scaled[x,y] <- matrix_1_vs_2.5[x,y]*scaling_factor
      }
    }
  }
  max_saved_blue <- paste("max saved blue = ", round(max(matrix_1_vs_2.5)))
  Palette_choice <- jBluGrnPalette
} else {
  matrix_scaled <- -matrix_scaled
  max_saved_green <- paste("max saved green = 0")
  Palette_choice <- jBluesPalette
}
heatmap(matrix_scaled, scale="none", Rowv=NA, Colv=NA, main=graph_title2, col = Palette_choice, xlab="v1_day", ylab="vac1")
mtext(max_saved_blue, 1, line = -12, adj = 1.6)
mtext(max_saved_green, 1, line = -11, adj = 1.8)
mtext(vacs, 1, line = -10, adj = 1.4)

#################################
# Different Reproductive Number #
#################################
parms <- parms_original
inits <- inits_original
vaccines <- 20000
bins2 <- 20
R0_timing <- data.frame(matrix(rep(0,bins2*2),nrow=bins2))
R0_trends <- data.frame(matrix(rep(0,bins2*4*5),nrow=bins2*4))
colnames(R0_timing) <- c("R0", "peak_time")
colnames(R0_trends) <- c("R0", "percent", "AB_1dose", "A_2dose", "AB_2dosehalf")

# calculate peak epidemic time
for (i in 1:bins2){
  parms["B0"] <- i/(.75*bins2)
  R0 <- parms["B0"]/parms["r"]
  R0_timing[i,1] <- R0
  parms["v1_count"] <- 0
  parms["v2_count"] <- 0
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  I_sum <- simulation.df[c("I_0")] + simulation.df[c("I_1")] + simulation.df[c("I_2")]
  R0_timing[i,2] <- which.max(I_sum[,1])
}
View(R0_timing)

times <- c(.25, .50, .75, 1.00)
R0_trends[,1] <- rep(R0_timing[,1],4)
for (j in 1:4){
  time <- times[j]
  R0_trends[((j-1)*bins2+1):(j*bins2),2] <- paste(time*100, "% of peak")
  # assume a fixed number of vaccines = enough to vaccinate everyone once
  for (i in 1:bins2){
    if (R0_timing[i,2] > 6){
      # an epidemic occurs -> record number of infections
      parms["v1_day"] <- round(R0_timing[i,2]*time,0)
      parms["v2_day"] <- parms["v1_day"] + 14
      # 1 dose AB
      parms["v1_count"] <- vaccines/2
      parms["v2_count"] <- 0
      simulation0 <- Run_Model_D(inits, dt, parms=parms)
      simulation0.df <- as.data.frame(simulation0)
      simulation0.df$D <- 1-apply(simulation0.df[,2:length(simulation0.df)], 1, sum)
      R0_trends[((j-1)*bins2+i),3] <- 2*sum(simulation0.df[nrow(simulation0.df), c("R_0", "R_1","R_2")])
      # 2 dose A
      parms["v1_count"] <- vaccines/2
      parms["v2_count"] <- vaccines/2
      simulation1a <- Run_Model_D(inits, dt, parms=parms)
      simulation1a.df <- as.data.frame(simulation1a)
      simulation1a.df$D <- 1-apply(simulation1a.df[,2:length(simulation1a.df)], 1, sum)
      R0_trends[((j-1)*bins2+i),4] <- sum(simulation1a.df[nrow(simulation1a.df), c("R_0", "R_1","R_2")]) 
      # 2 dose half AB
      parms["v1_count"] <- vaccines/4
      parms["v2_count"] <- vaccines/4
      simulation2 <- Run_Model_D(inits, dt, parms=parms)
      simulation2.df <- as.data.frame(simulation2)
      simulation2.df$D <- 1-apply(simulation2.df[,2:length(simulation2.df)], 1, sum)
      R0_trends[((j-1)*bins2+i),5] <- 2*sum(simulation2.df[nrow(simulation2.df), c("R_0", "R_1","R_2")])
    }
  }
}
# 0 dose (same for all v1_days and v2_days)
parms["v1_count"] <- 0
parms["v2_count"] <- 0
simulation1b <- Run_Model_D(inits, dt, parms=parms)
simulation1b.df <- as.data.frame(simulation1b)
simulation1b.df$D <- 1-apply(simulation1b.df[,2:length(simulation1b.df)], 1, sum)

for (k in 1:nrow(R0_trends)){
  if (k%%nrow(R0_timing) != 0){
    m <- k%%nrow(R0_timing)
  } else{
    m <- bins2
  }
    if (R0_timing[m,2] > 6){
      R0_trends[k,4] <- R0_trends[k,4] + sum(simulation1b.df[nrow(simulation1b.df), c("R_0", "R_1","R_2")])
    }
}
View(R0_trends)
ggplot(data=R0_trends, aes(x=R0)) + 
  geom_line(aes(x=R0, y=AB_1dose), col="blue") + 
  geom_line(aes(x=R0, y=A_2dose), col="red") + 
  geom_line(aes(x=R0, y=AB_2dosehalf), col="green") + 
  facet_wrap( ~ percent, ncol=2) + 
  ggtitle("Epidemic Size According to R0")

########################################
# Incremental Benefit of a Second Dose #
########################################
# vary the percent coverage of a population for first and second doses
bins <- 11
parms <- parms_original
inits <- inits_original
add_lives_saved <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
percentcoverage_round1 <- c(rep(NA,bins))
percentcoverage_round2 <- c(rep(NA,bins))
for (i in 1:bins){
  percentcoverage_round1[i] <- (i-1)*10
  parms["v1_count"] <- parms["ntot"]*percentcoverage_round1[i]/100
  for (j in 1:i){
    percentcoverage_round2[j] <- (j-1)*10
    parms["v2_count"] <- parms["ntot"]*percentcoverage_round2[j]/100
    simulation <- Run_Model_D(inits, dt, parms=parms)
    simulation.df <- as.data.frame(simulation)
    simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
    add_lives_saved[i,j] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])
  }
}
rownames(add_lives_saved) <- percentcoverage_round1
colnames(add_lives_saved) <- percentcoverage_round2

# # compare to max lives lost (i.e. no coverage)
# max_lives_lost <- max(add_lives_saved)
# add_lives_saved_ml <- add_lives_saved
# for (x in 1:bins){
#   for (y in 1:x){
#     add_lives_saved_ml[x,y] <- max_lives_lost - add_lives_saved_ml[x,y]
#   }
# }
## flip to epidemic size NO RED
# als_ml <- t(as.matrix(add_lives_saved_ml))
als <- t(as.matrix(add_lives_saved))
graph_title2 <- "epidemic size"
Palette_choice <- jBluesPalette
heatmap(als, scale="none", Rowv=NA, Colv=NA, main=graph_title2, col = Palette_choice, xlab="percent coverage of first dose", ylab="percent coverage of second dose")

# compare to lives lost for each percent coverage of first dose
add_lives_saved_pc1 <- add_lives_saved
for (x in 1:bins){
  for (y in 1:x){
    add_lives_saved_pc1[x,y] <- add_lives_saved[x,1] - add_lives_saved_pc1[x,y]
  }
}
als_pc1 <- t(as.matrix(add_lives_saved_pc1))
graph_title2 <- "Additional lives saved compared to 0 second dose"
heatmap(als_pc1, scale="none", Rowv=NA, Colv=NA, main=graph_title2, col = Palette_choice, xlab="percent coverage of first dose", ylab="percent coverage of second dose")

# above, but per second dose administered
add_lives_saved_per_vac <- add_lives_saved_pc1
for (x in 1:bins){
  for (y in 1:x){
    add_lives_saved_per_vac[x,y] <- add_lives_saved_per_vac[x,y]/(parms["ntot"]*(y-1)*100)
  }
}
als_pc1_vac <- t(as.matrix(add_lives_saved_per_vac))
graph_title2 <- "ALS per 2nd dose compared to 0 second dose"
heatmap(als_pc1_vac, scale="none", Rowv=NA, Colv=NA, main=graph_title2, col = Palette_choice, xlab="percent coverage of first dose", ylab="percent coverage of second dose")

##
## vary the timing of the first dose and its effectiveness
bins <- 21
parms <- parms_original
parms["v1_count"] <- .9*parms["ntot"]
parms["v2_count"] <- .8*parms["ntot"]
inits <- inits_original
add_lives_saved <- data.frame(matrix(rep(NA, bins*bins), nrow=bins))
timing_1 <- c(rep(NA,bins))
vac_1 <- c(rep(NA,bins))
for (i in 1:bins){
  timing_1[i] <- i*10
  parms["v1_day"] <- timing_1[i]
  parms["v2_day"] <- parms["v1_day"] + 14
  for (j in 1:bins){
    vac_1[j] <- (j-1)/(bins-1) * .5
    parms["vac1"] <- vac_1[j]
    simulation <- Run_Model_D(inits, dt, parms=parms)
    simulation.df <- as.data.frame(simulation)
    simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
    add_lives_saved[i,j] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])
  }
}
rownames(add_lives_saved) <- timing_1
colnames(add_lives_saved) <- vac_1

max_saved <- 14571
als <- t(as.matrix(max_saved - add_lives_saved))
graph_title2 <- "Cases averted"
Palette_choice <- jGreensPalette
heatmap(als, scale="none", Rowv=NA, Colv=NA, main=graph_title2, col = Palette_choice, xlab="timing of first dose", ylab="vac1_efficacy")
View(als)

# compare to lives lost for each percent coverage of first dose
add_lives_saved_pc1 <- add_lives_saved
for (x in 1:bins){
  for (y in 1:x){
    add_lives_saved_pc1[x,y] <- add_lives_saved[x,1] - add_lives_saved_pc1[x,y]
  }
}
als_pc1 <- t(as.matrix(add_lives_saved_pc1))
graph_title2 <- "Additional lives saved compared to 0 second dose"
heatmap(als_pc1, scale="none", Rowv=NA, Colv=NA, main=graph_title2, col = Palette_choice, xlab="percent coverage of first dose", ylab="percent coverage of second dose")

# above, but per second dose administered
add_lives_saved_per_vac <- add_lives_saved_pc1
for (x in 1:bins){
  for (y in 1:x){
    add_lives_saved_per_vac[x,y] <- add_lives_saved_per_vac[x,y]/(parms["ntot"]*(y-1)*100)
  }
}
als_pc1_vac <- t(as.matrix(add_lives_saved_per_vac))
graph_title2 <- "ALS per 2nd dose compared to 0 second dose"
heatmap(als_pc1_vac, scale="none", Rowv=NA, Colv=NA, main=graph_title2, col = Palette_choice, xlab="percent coverage of first dose", ylab="percent coverage of second dose")

#############################
# Effectives and Efficiency #
#############################

parms <- parms_original
inits <- inits_original
coverage_vac1 <- .90
coverage_vac2 <- .80
end <- 150
sensitivity <- data.frame(matrix(rep(0,end*4),nrow=end))
ep_size <- data.frame(matrix(rep(0,end*8),nrow=end))
ep_size[,1] <- 1:end
colnames(ep_size) <- c("v1_day", "nodose", "dose1", "dose2", "delta_21", "delta_21rel", "delta_01", "delta_02")
for (i in 3:end){
  parms["v1_day"] <- i
  parms["v2_day"] <- parms["v1_day"] + 14
  # no dose
  parms["v1_count"] <- 0
  parms["v2_count"] <- 0
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  ep_size[i,2] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])
  # 1st dose only
  parms["v1_count"] <- coverage_vac1*parms["ntot"]
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  ep_size[i,3] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])
  # 2nd dose
  parms["v2_count"] <- coverage_vac2*parms["ntot"]
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  ep_size[i,4] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")]) 
}
# ggplot(data=ep_size, aes(x=v1_day))+
#   geom_line(aes(x=v1_day, y=dose1), col="red")+
#   geom_line(aes(x=v1_day, y=dose2), col="green")+
#   xlab("v1_day") + ylab("final epidemic size (# people)")+
#   ggtitle("Absolute Effectiveness")

ep_size[,5] <- ep_size[,3] - ep_size[,4]
ggplot(data=ep_size, aes(x=v1_day))+
  geom_line(aes(x=v1_day, y=delta_21), col="black")+
  xlab("v1_day") + ylab("cases averted by using 2nd dose")+
  ggtitle("Relative Effectiveness")

ep_size[,6] <- ep_size[,5]/ep_size[,3]*100
ggplot(data=ep_size, aes(x=v1_day))+
  geom_line(aes(x=v1_day, y=delta_21rel), col="black")+
  xlab("v1_day") + ylab("100*cases averted using 2nd dose/final ep size of 1st dose")+
  ggtitle("Relative Effectiveness (as a percent, compared to # averted by 1st dose)")

# cases averted compared to no dose per thousand doses (div 9 and 8)
ep_size[,7] <- (ep_size[,2] - ep_size[,3])/(coverage_vac1*1000)
ep_size[,8] <- (ep_size[,2] - ep_size[,4])/(coverage_vac1*1000+coverage_vac2*1000)

ep_efficiency <- ep_size[,6:8]
ep_efficiency[,1] <- ep_size[,1]
colnames(ep_efficiency) <- c("v1_day", "dose1", "dose2")
ep_eff_m <- melt(ep_efficiency, id.var="v1_day")
ggplot(data=ep_eff_m, aes(x=v1_day, y=value))+
  geom_line(aes(colour=variable, group=variable))+
  xlab("v1_day") + ylab("cases averted per thousand doses")+
  ggtitle("Absolute Efficiency")

rel_eff <- ep_efficiency[,2] - ep_efficiency[,3]
rel <- ep_size[,1:2]
rel[,2] <- rel_eff
colnames(rel) <- c("v1_day", "diff")
ggplot(data=rel, aes(x=v1_day))+
  geom_line(aes(x=v1_day, y=diff), col="black")+
  xlab("v1_day") + ylab("cases averted per thousand doses")+
  ggtitle("Relative Efficiency")


# sensitivity analysis
ep2_size <- ep_size[,1:4]
ep3_size <- ep_size[,1:4]

parms["vac1"] <- .3
for (i in 3:end){
  parms["v1_day"] <- i
  parms["v2_day"] <- parms["v1_day"] + 14
  # no dose
  parms["v1_count"] <- 0
  parms["v2_count"] <- 0
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  ep2_size[i,2] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])
  # 1st dose only
  parms["v1_count"] <- coverage_vac1*parms["ntot"]
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  ep2_size[i,3] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])
  # 2nd dose
  parms["v2_count"] <- coverage_vac2*parms["ntot"]
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  ep2_size[i,4] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")]) 
}

parms["vac1"] <- .7
for (i in 3:end){
  parms["v1_day"] <- i
  parms["v2_day"] <- parms["v1_day"] + 14
  # no dose
  parms["v1_count"] <- 0
  parms["v2_count"] <- 0
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  ep3_size[i,2] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])
  # 1st dose only
  parms["v1_count"] <- coverage_vac1*parms["ntot"]
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  ep3_size[i,3] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])
  # 2nd dose
  parms["v2_count"] <- coverage_vac2*parms["ntot"]
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  ep3_size[i,4] <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")]) 
}

sensitivity[,1] <- 1:end
sensitivity[,2] <- ep_size[,3] - ep_size[,4]
sensitivity[,3] <- ep2_size[,3] - ep2_size[,4]
sensitivity[,4] <- ep3_size[,3] - ep3_size[,4]
colnames(sensitivity) <- c("v1_day", "five", "three", "seven")

sensitivity[1:2,3] <- NA
ggplot(data=sensitivity, aes(x=v1_day))+
  geom_line(aes(x=v1_day, y=five), col="orange")+
  geom_line(aes(x=v1_day, y=three), col="purple")+
  geom_line(aes(x=v1_day, y=seven), col="dark green")+
  xlab("v1_day") + ylab("cases averted by using 2nd dose")+
  ggtitle("Relative Effectiveness")

################
# Effect of R0 #
################

parms <- parms_original
inits <- inits_original
simulation <- Run_Model_D(inits, dt, parms=parms)
simulation.df <- as.data.frame(simulation)
simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)

############################
# Model Performance Graphs #
############################

# check model for varying infectiousness
parms["v1_count"] <- 0
parms["v2_count"] <- 0
parms["B0"] <- .4
simulation1 <- Run_Model_D(inits, dt, parms=parms)
simulation1.df <- as.data.frame(simulation1)
simulation1.df$D <- 1-apply(simulation1.df[,2:length(simulation1.df)], 1, sum)

compare <- matrix(rep(0,3*nrow(simulation1)),nrow(simulation1))
compare[,1] <- simulation1.df[,1]
compare[,2] <- simulation1.df[,11]

parms["v1_count"] <- 0
parms["v2_count"] <- 0
parms["B0"] <- .8
simulation2 <- Run_Model_D(inits, dt, parms=parms)
simulation2.df <- as.data.frame(simulation2)
simulation2.df$D <- 1-apply(simulation2.df[,2:length(simulation2.df)], 1, sum)
compare[,3] <- simulation2.df[,11]
colnames(compare) <- c("time", "I_0_low", "I_0_high")
compare.df <- as.data.frame(compare)
ggplot(data=compare.df, aes(x=time)) +
  geom_line(aes(x=time, y=I_0_low), col="blue") +
  geom_line(aes(x=time, y=I_0_high), col="green") +
  theme_bw() + ylab("Count") + ggtitle("SEIR")

# vaccine introduction
parms <- parms_original
parms["B0"] <- .64
parms["v1_count"] <- .2*parms["ntot"]
parms["v2_count"] <- .2*parms["ntot"]
parms["v1_day"] <- 30
parms["v2_day"] <- 70
simulation <- Run_Model_D(inits, dt, parms=parms)
simulation.df <- as.data.frame(simulation)
simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
Plot_Compartments(simulation.df)

ggplot(data=simulation.df, aes(x=time)) +
  geom_line(aes(x=time, y=S_0), col="blue") +
  geom_line(aes(x=time, y=S_1), col="green") +
  geom_line(aes(x=time, y=S_2), col="black") +
  theme_bw() + ylab("Count") + ggtitle("SEIR")


ggplot(data=simulation.df, aes(x=time)) +
  geom_line(aes(x=time, y=S_0), col="blue", linetype="dotted") +
  geom_line(aes(x=time, y=S_1), col="blue", linetype="dashed") +
  geom_line(aes(x=time, y=S_2), col="blue", linetype="solid") +
  theme_bw() + ylab("Count") + ggtitle("SEIR")

######################
# Epidemic Size, 3/7 #
######################
parms <- parms_original
t_final <- 215
inits <- c(S_0=24999,
           S_1=0,
           S_2=0,
           E_0=0,
           E_1=0,
           E_2=0,
           I_0=1,
           I_1=0,
           I_2=0,
           R_0=0,
           R_1=0,
           R_2=0,
           N  =25000)
start <- 6
delay <- 14 #Note, Corey changed this to 14 on Dec 14
end <- t_final - delay - 1
time_steps <- end - start + 1
epidemic_size <- data.frame(matrix(rep(0,time_steps*4),nrow=time_steps))
epidemic_size[,1] <- start:end
colnames(epidemic_size) <- c("time", "dose0", "dose1", "dose2")
# assume a fixed number of vaccines = enough to vaccinate everyone once
for (i in start:end){
  parms["v1_day"] <- i
  parms["v2_day"] <- parms["v1_day"] + delay
  # 1 dose 
  parms["v1_count"] <- .9*inits["N"]
  parms["v2_count"] <- 0
  simulation0 <- Run_Model_D(inits, dt, parms=parms)
  simulation0.df <- as.data.frame(simulation0)
  simulation0.df$D <- 1-apply(simulation0.df[,2:length(simulation0.df)], 1, sum)
  epidemic_size[i,3] <- sum(simulation0.df[nrow(simulation0.df), c("R_0", "R_1","R_2")])
  # 2 dose
  parms["v1_count"] <- .9*inits["N"]
  parms["v2_count"] <- .8*inits["N"]
  simulation1a <- Run_Model_D(inits, dt, parms=parms)
  simulation1a.df <- as.data.frame(simulation1a)
  simulation1a.df$D <- 1-apply(simulation1a.df[,2:length(simulation1a.df)], 1, sum)
  epidemic_size[i,4] <- sum(simulation1a.df[nrow(simulation1a.df), c("R_0", "R_1","R_2")]) 
  }
# 0 dose (same for all v1_days and v2_days)
parms["v1_count"] <- 0
parms["v2_count"] <- 0
simulation1b <- Run_Model_D(inits, dt, parms=parms)
simulation1b.df <- as.data.frame(simulation1b)
simulation1b.df$D <- 1-apply(simulation1b.df[,2:length(simulation1b.df)], 1, sum)
epidemic_size[,2] <- c(rep(sum(simulation1b.df[nrow(simulation1b.df), c("R_0", "R_1","R_2")]), end))
epidemic_size <- epidemic_size[1:end,]
ggplot(data=epidemic_size, aes(x=time)) + 
  geom_line(aes(x=time, y=dose0), col="blue") + 
  geom_line(aes(x=time, y=dose1), col="red") + 
  geom_line(aes(x=time, y=dose2), col="green") + 
  theme_bw() + xlab("Time Until First Vaccine") + ylab("Cumulative Epidemic Size") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Epidemic Size According to Campaign")

# plot epidemic with vaccination after 13th week
parms <- parms_original
parms["v1_day"] <- 200
parms["v2_day"] <- parms["v1_day"] + 14
parms["v1_count"] <- .9*parms["ntot"]
parms["v2_count"] <- .8*parms["ntot"]
simulation <- Run_Model_D(inits, dt, parms=parms)
simulation.df <- as.data.frame(simulation)
simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
simulation.df <- simulation.df[1:200,]
Plot_SEIR(simulation.df)
before <- sum(simulation.df[parms["v1_day"]-1, c("I_0", "I_1", "I_2", "R_0", "R_1","R_2")])
before1 <- sum(simulation.df[parms["v1_day"]-1, c("R_0", "R_1","R_2")])
after <- sum(simulation.df[nrow(simulation.df), c("R_0", "R_1","R_2")])
diff <- after - before

eff <- epidemic_size[,3] - epidemic_size[,4]

# calculate change
change <- matrix(rep(0,2*(nrow(epidemic_size)-1)),(nrow(epidemic_size)-1))
for (i in 1:(nrow(epidemic_size)-1)){
  change[i,1] <- epidemic_size[i+1,3] - epidemic_size[i,3]
  change[i,2] <- epidemic_size[i+1,4] - epidemic_size[i,4]
}

##########################
# Epidemic Duration, 3/7 #
##########################
parms <- parms_original
parms["vac1"] <- .65
t_final <- 200
inits <- c(S_0=24999,
           S_1=0,
           S_2=0,
           E_0=0,
           E_1=0,
           E_2=0,
           I_0=1,
           I_1=0,
           I_2=0,
           R_0=0,
           R_1=0,
           R_2=0,
           N  =25000)
start <- 6
delay <- 14 #Note, Corey changed this to 14 on Dec 14
end <- t_final - delay - 1
time_steps <- end - start + 1
epidemic_duration <- data.frame(matrix(rep(0,time_steps*4),nrow=time_steps))
epidemic_duration[,1] <- start:end
colnames(epidemic_duration) <- c("time", "dose0", "dose1", "dose2")
for (i in start:end){
  parms["v1_day"] <- i
  parms["v2_day"] <- parms["v1_day"] + delay
  # 1 dose 
  parms["v1_count"] <- .9*inits["N"]
  parms["v2_count"] <- 0
  simulation0 <- Run_Model_D(inits, dt, parms=parms)
  simulation0.df <- as.data.frame(simulation0)
  simulation0.df$D <- 1-apply(simulation0.df[,2:length(simulation0.df)], 1, sum)
  count <- 1
  while(.99 < sum(simulation0.df[count, c("E_0", "E_1","E_2", "I_0", "I_1","I_2")])) {
    count <- count + 1
  }
  epidemic_duration[i,3] <- simulation0.df[count,1]
  # 2 dose
  parms["v1_count"] <- .9*inits["N"]
  parms["v2_count"] <- .8*inits["N"]
  simulation1a <- Run_Model_D(inits, dt, parms=parms)
  simulation1a.df <- as.data.frame(simulation1a)
  simulation1a.df$D <- 1-apply(simulation1a.df[,2:length(simulation1a.df)], 1, sum)
  count <- 1 #Corey added this on March 23
  while(0.99 < sum(simulation1a.df[count, c("E_0", "E_1","E_2", "I_0", "I_1","I_2")])) {
    count <- count + 1
  }
  epidemic_duration[i,4] <- simulation1a.df[count,1]
}
# 0 dose (same for all v1_days and v2_days)
parms["v1_count"] <- 0
parms["v2_count"] <- 0
simulation1b <- Run_Model_D(inits, dt, parms=parms)
simulation1b.df <- as.data.frame(simulation1b)
simulation1b.df$D <- 1-apply(simulation1b.df[,2:length(simulation1b.df)], 1, sum)
count <- 1 #Corey added this on March 23
while(0.99 < sum(simulation1b.df[count, c("E_0", "E_1","E_2", "I_0", "I_1","I_2")])) {
  count <- count + 1
}
epidemic_duration[,2] <- c(rep(simulation1b.df[count,1],end))
epidemic_duration <- epidemic_duration[6:180,]

ggplot(data=epidemic_duration, aes(x=time)) + 
  geom_line(aes(x=time, y=dose0), col="blue") + 
  geom_line(aes(x=time, y=dose1), col="red") + 
  geom_line(aes(x=time, y=dose2), col="green") + 
  theme_bw() + xlab("Time Until First Vaccine (days)") + ylab("Epidemic Duration (days)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Epidemic Duration According to Campaign")

Plot_SEIR(simulation0.df)

#######################
# Compare delay times #
#######################
bins <- 12
parms <- parms_original
inits <- inits_original
parms["v1_count"] <- parms["ntot"]*.9
parms["v2_count"] <- parms["ntot"]*.8
heat_map_timing <- data.frame(matrix(rep(NA, bins*bins), nrow=bins))
rownames(heat_map_timing) <- as.character((1:bins)*bins)

for (x in 1:bins){
  parms["v1_day"] <- 10*x
  parms["v2_day"] <- 350
  rownames(heat_map_timing)[x] <- parms["v1_day"]
  simulation <- Run_Model_D(inits, dt, parms=parms)
  simulation.df <- as.data.frame(simulation)
  simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
  cases <- sum(simulation.df[nrow(simulation.df),c("R_0", "R_1", "R_2")])
  for (y in 1:bins){
    parms["v2_day"] <- parms["v1_day"] + 2*y
    colnames(heat_map_timing)[y] <- 2*y
    simulation0 <- Run_Model_D(inits, dt, parms=parms)
    simulation0.df <- as.data.frame(simulation0)
    simulation0.df$D <- 1-apply(simulation0.df[,2:length(simulation0.df)], 1, sum)
    cases2 <- sum(simulation0.df[nrow(simulation0.df), c("R_0", "R_1", "R_2")])
    heat_map_timing[x,y] <- cases - cases2
  }
}
hm_timing <- as.matrix(heat_map_timing)
# create a heatmap 
heatmap(hm_timing, Rowv = NA, Colv = NA, scale = "none", main = "v1_day vs. v2_day", col = jGreensPalette)

# sanity check that heatmap is accurate by examining 2 week delay only
heat_map14 <- heat_map_timing
heat_map14[,1:6] <- NA
heat_map14[,8:12] <- NA
hm14 <- as.matrix(heat_map14)
heatmap(hm14, Rowv = NA, Colv = NA, scale = "none", main = "14 day delay", col = jGreensPalette)

##################################
# Number of Doses Allocated 3/25 #
##################################
bins <- 13
parms <- parms_original # reset parameters
parms["v1_day"] <- 50
parms["v2_day"] <- parms["v1_day"] + 14
heat_map_reac <- data.frame(matrix(rep(0, bins*bins), nrow=bins))
rownames(heat_map_reac) <- as.character((1:bins)/100) #just a work-around to initialize
for (x in 1:bins){
  parms["v1_count"] = 2000*x - 1000
  rownames(heat_map_reac)[x] <- as.character(parms["v1_count"])
  for (y in 1:bins){
    parms["v2_count"] = 2000*y - 1000
    colnames(heat_map_reac)[y] <- as.character(parms["v2_count"])
    simulation <- Run_Model_D(inits, dt, parms=parms)
    simulation.df <- as.data.frame(simulation)
    simulation.df$D <- 1-apply(simulation.df[,2:length(simulation.df)], 1, sum)
    heat_map_reac[x,y] <- sum(simulation.df[nrow(simulation.df),c("R_0", "R_1", "R_2")])
  }
}

# create a heatmap 
matrix_reac <- as.matrix(t(heat_map_reac))
heatmap(matrix_reac, Rowv = NA, Colv = NA, scale = "none", main = "v1_count vs. v2_count", xlab = "v1_count", ylab = "v2_count", col = jBluesPalette)
