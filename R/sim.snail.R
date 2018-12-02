
sim.snail<-function(MaxX=53, MinX=0, MaxY=7.5, MinY=0, Pstop=0.05, Steps=28800,
                     Vels=Vels, b=15, m=0.3018, Q10=2, Tref=22,
                     mua=mua, rhoa=rhoa, rhostop=0.5*rhoa,
                     rhowall=0.9, export="meanlasthour"){

  #library(circular)
  library(CircStats)
  
  # Starting Conditions
  StartX<-(MaxX-MinX)/2
  StartY<-(MaxY-MinY)/2
  
  # Draw a random starting direction
  # ThetaStart<-rwrpnorm(1, mu=mut, rho=rhot)
  ThetaStart<-runif(1, min=0, max=2*pi)
  
  # Create a multimodal distribution to deal with what angle to take when bouncing
  # off a wall:
  
  horizangles<-c(rwrpnorm(10000, mu=0, rho=rhowall), rwrpnorm(10000, mu=pi, rho=rhowall))
  vertangles<-c(rwrpnorm(10000, mu=pi/2, rho=rhowall), rwrpnorm(10000, mu=pi*3/2, rho=rhowall))
  
  X<-rep(NA, Steps); Y<-rep(NA, Steps)
  X[1]<-StartX; Y[1]<-StartY
  # Z<-rep(NA, Steps); Z[1]<-0
  t<-rep(NA, Steps); t[1]<-ThetaStart
  a<-rep(NA, Steps); a[1]<-0
  Temperature<-rep(NA, Steps); Temperature[1]<-b+m*StartX
  rhos<-c(rhostop, rhoa)
  pmove<-rbinom(Steps, 1, prob=1-Pstop) # prepopulate
  rhoalpha<-rep(rhoa,Steps)   # create a vector where rhoalpha = rhoa
  rhoalpha[pmove==0]<-rhostop # then set each index where pmove=0 to rhostop
  # returns a working rho value depending on the indexed pmove
  # then calculate na, using this working rho.  I'd like to minimise
  # using the if statement below
  
  foo<-function(mua, rhoalpha) return(rwrpnorm(1, mua, rho=rhoalpha))
  na<-sapply(X=rhoalpha, foo, mua=mua) # prepopulated
  
  # plot(StartX, StartY, xlim=c(MinX, MaxX), ylim=c(MinY, MaxY))
  
  for(i in 2:Steps){
    
    dz<-sample(Vels, 1) # draw a random sample from the Velocity observations
    
    dz<-dz * 10^(log10(Q10)*(Temperature[i-1]-Tref)/10) 
    # Q10 correct the movement rate using Tref=22C as reference temperature
    
    # rhoalpha<-rhos[pmove[i]+1]  # vectorised and slower approach
    
    # na<-rwrpnorm(1, mua, rhoalpha[i]) # vectoried and slower approach
    
    # if pmove=0, then calculate turning angle to draw from a more variable distribution
    # where rho is rhostop as determined in function call.  Probably useful to 
    # make rhostop = 0.5*rhoa (wider variance)
    
    # # if Pmove =1, then calculate turning angle, na, normally
    # if(pmove[i]==1) {
    #   na<-rwrpnorm(1, mua, rhoa)
    # } else {
    #   na<-rwrpnorm(1, mua, rho=rhostop) # if(pmove[i]==0)  
    # } 
    
    nt<-t[i-1] + na[i]
    
    #dx<-cos(nt)*dz # i.e. adjacent<-cos(Theta)*hypotenuse (where dz=hypotenuse)
    #dy<-sin(nt)*dz # i.e. opposite<-sin(Theta)*hypotenuse
    
    X[i] <- X[i-1] + cos(nt)*dz*pmove[i]
    Y[i] <- Y[i-1] + sin(nt)*dz*pmove[i]
    
    # If outside the boundaries, recalculate the angle with less
    # stringent autocorrelation:
    
    while(X[i] < MinX | X[i] > MaxX | Y[i] < MinY | Y[i] > MaxY) {
      # if the walk falls outside the boundaries, redraw a new theta
      # from a circular normal distribution with mean similar to 
      # observed snail movements, but with a rho that has high 
      # variability
      
      if(X[i] < MinX | X[i] > MaxX) nt<-sample(vertangles, 1) # draw from the bimodal vertangles
      if(Y[i] < MinY | Y[i] > MaxY) nt<-sample(horizangles, 1) # draw from the bimodal horizangles
      
      #dx<-cos(nt)*dz # i.e. adjacent<-cos(Theta)*hypotenuse (where dz=hypotenuse)
      #dy<-sin(nt)*dz # i.e. opposite<-sin(Theta)*hypotenuse
      
      X[i] <- X[i-1] + cos(nt)*dz*pmove[i]
      Y[i] <- Y[i-1] + sin(nt)*dz*pmove[i]
      
      nt<-atan2(Y[i]-Y[i-1], X[i]-X[i-1])
    } 
    
    
    Temperature[i]<-b+m*X[i]
    
    t[i]<-nt
    # a[i]<-t[i]-t[i-1]
    # Z[i]<-dz
    
  }
  
  if(export=="meanlasthour") output<-mean(Temperature[(Steps-3600):Steps])
  if(export=="path") output<-data.frame(X,Y,t,pmove,Temperature)
  return(output)
}


rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}
CI<-function(x) 1.96*sd(x)/sqrt(length(x))
lower25<-function(x) quantile(x, 0.25)
upper75<-function(x) quantile(x, 0.75)
circmean<-function(x) circ.mean(na.omit(x))
circrho<-function(x) circ.disp(na.omit(x))$rbar
