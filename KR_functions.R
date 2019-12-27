########## set initial conditions, fitness costs, and free parameters
set_init <- function(experiment,genotype_names){
  init <- rep(0,length(genotype_names))
  eq <- 25
  init[grepl("kkrr",genotype_names)] <- eq
  ## Set fitness values - parameters to be estimated will be changed to 0.9
  w <- rep(1,length(genotypes))
  w[grepl("KKrr",genotype_names)] <- 0 #no rescue always lethal
  w[grepl("Kkrr",genotype_names)] <- 0
  w[grepl("KKRR",genotype_names)] <- 1
  w[grepl("KkRR",genotype_names)] <- 1
  w[grepl("KKRr",genotype_names)] <- 1
  w[grepl("KkRr",genotype_names)] <- 1
  w[grepl("kkRr",genotype_names)] <- 1
  w[grepl("kkRR",genotype_names)] <- 0.87
  
  names(w) <- genotype_names
  #initial conditions differ by experiment. pindex 
  if (experiment == "DG4") {
    init[grepl("KKRR",genotype_names)] <- 2*eq
    pindex<-c(grep("KKRR",genotype_names),grep("KKRr",genotype_names),
              grep("KkRR",genotype_names),grep("KkRr",genotype_names))
    tend <- 9
  } else if (experiment == "HG4" | experiment=="HG4hid"){
    init[grepl("KkRR",genotype_names)] <- 2*eq
    w[grepl("KKRR",genotype_names)] <- 0
    w[grepl("KKRr",genotype_names)] <- 0
    pindex<-c(grep("KkRR",genotype_names),grep("KkRr",genotype_names),
              grep("kkRR",genotype_names),grep("kkRr",genotype_names))
    tend <- 6
  } else if (experiment == "G80"){
    init[grepl("kkRR",genotype_names)] <- 2*eq
    pindex<-c(grep("kkRR",genotype_names),grep("kkRr",genotype_names))
    tend <- 6
  }
  list(init=init,w=w,pindex=pindex,tendin=tend,experiment=experiment,genotype_names=genotype_names)
}

########## deterministic simulation
simfunc <- function(input) {
  Parray <- input$Parray
  ngenos <- dim(Parray)[1] #number of genotypes = 9
  
  tend <- input$tendin
  A <- matrix(nrow=tend,ncol=ngenos,dimnames = list(NULL,input$genotype_names))
  B <- A
  
  #initial conditions from input
  A[1,] <- input$init/sum(input$init)
  
  alloff <- array(dim=c(tend,dim(Parray)))
  for(t in 2:(tend)){
    for(i in 1:(ngenos)){
      for(m in 1:(ngenos)){
        alloff[t,i,m,] <- input$w[i]*A[t-1,m]*(1)*Parray[m,,i]*A[t-1,]
      }
      B[t,i] <- sum(alloff[t,i,,])
    }
    A[t,] <- B[t,]/sum(B[t,])
  }
  simdat <- as.data.frame(A)/rowSums(A)
  simdat <- as.data.frame(A)
  simdat <- simdat %>% subset(select=-c(KKrr,Kkrr))
  allsimdat <- simdat %>% within({
    time <- 1:tend
    replicate <- "0"
    experiment <- input$experiment
    simtype <- "sim"
  })
  allsimdat
}

########## stochastic simulation
# equivalent calculations for the deterministic model are included as comments
simfunc_stoch <- function(input,rep=1){
  Parray <- input$Parray
  
  #automatically make ordering of genotypes from above, and assign fitness cost vectors accordingly
  ngenos <- dim(Parray)[1]
  tend <- input$tendin
  females <- matrix(nrow=tend,ncol=ngenos,dimnames = list(NULL,input$genotype_names))
  females[1,] <- 0
  males <- females
  B <- females
  females[1,] <- input$init
  males[1,] <- input$init
  total <- females+males

  lambda_base <- 6
  lambda <- rep(lambda_base,ngenos) #assume all genotypes produce same mean number of offspring
  
  alloff <- array(dim=c(tend,ngenos,ngenos,ngenos))
  
  classes <- array(dim=c(tend,ngenos,ngenos))
  classeggs <- classes
  for(t in 2:(tend)){
    if (sum(males[t-1,])==0 & sum(females[t-1,])==0){ #check to make sure we're mating adults that exist
      alloff[t,,,] = 0
    } else{
      for(m in 1:ngenos){
        matingfemales <- females[t-1,m] #rbinom(1,females[t-1,m],1-whiteeye[m]*dwh) #randomly draw to check surviving females
        classes[t,m,] <- rmultinom(1,size=matingfemales,prob = males[t-1,]/sum(males[t-1,]))
        # classes[t,m,] <- round(matingfemales* males[t-1,]/sum(males[t-1,]))
        
        classeggs[t,m,] <- rpois(ngenos,classes[t,m,]*lambda[m]) #eggs per male mates
        # classeggs[t,m,] <- classes[t,m,]*lambda[m] #eggs per male mates
        
        for(n in 1:ngenos){
          alloff[t,,m,n] <- rmultinom(1,classeggs[t,m,n],Parray[m,n,])
          # alloff[t,,m,n] <- round(classeggs[t,m,n]*Parray[m,n,])
        }
      }
    }
    
    for(i in 1:ngenos){
      B[t,i] <- rbinom(1,sum(alloff[t,i,,]),input$w[i]) #tally and draw larvae that emerge (fitness costs)
      # B[t,i] <- round(sum(alloff[t,i,,])*input$w[i])
      females[t,i] <- rbinom(1,B[t,i],0.5) #[t,]*100/sum(B[t,])
      # females[t,i] <- round(B[t,i]*0.5) 
      males[t,i] <- B[t,i]-females[t,i]
    }
    adultcounts <- sum(females[t,]+males[t,])
    # adults <- females[t,]+males[t,]
    
    if(adultcounts>=(sum(input$init)*2)){
      # mysamples <- rmvhyper(1,c(females[t,],males[t,]),sum(input$initF)*2) #randomsampling
      mysamples <- round(sum(input$init)*2 / adultcounts * c(females[t,],males[t,])) #nonrandom
    } else {
      # mysamples <- rmvhyper(1,c(females[t,],males[t,]),adultcounts) #randomsampling
      mysamples <- c(females[t,],males[t,]) #nonrandom
    }
    females[t,] <- mysamples[1:ngenos]
    males[t,] <- mysamples[(ngenos+1):(2*ngenos)]
    total[t,] <- (females[t,]+males[t,])
  }
  
  simdat <- as.data.frame(total) #total instead of females
  simdat <- as.data.frame(total)/rowSums(total) #genotypes fractions instead of counts
  simdat <- simdat %>% subset(select=-c(KKrr,Kkrr)) #remove genotypes that are always zero
  allsimdat <- simdat %>% within({
    # sumJ <- rowSums(B)
    # sumF <- rowSums(females)
    # sumM <- rowSums(males)
    experiment <- input$experiment
    replicate <- as.character(rep)
    time <- 1:tend
    simtype <- "stochsim"
  })
  
  allsimdat
}

########## least squares function
KR_costfunc <- function(simdat,dat){
  x1 <- simdat[,unique(c(grep("K",names(simdat)),grep("k",names(simdat))))]
  x2 <- dat[,unique(c(grep("K",names(dat)),grep("k",names(dat))))]
  diffs <- (x1-x2[,names(x1)])
  sum(diffs^2)
}

########## minimizer function
KR_val_func <- function(p,input,dat){
  p <- transfrom(p)
  input$w[input$pindex] <- p
  simdat <- simfunc(input)
  simdat <- simdat[rep(seq(1,nrow(simdat)),times= length(unique(dat$replicate))),] #repeat data for each replicate
  cost <- KR_costfunc(simdat,dat) #calculate least squares
  cost
}

########## transform parameters for unconstrained optimization from [0,1]
transto <- function(p) (log(p/(1-p)))
transfrom <- function(p) 1/(1+exp(-p))

########## repeat this estimation function for parametric bootstrapping
bootfit <- function(input){  
  p <- rep(0.9,length(input$pindex))
  bootdat <- do.call(rbind, lapply(as.list(seq(1,5)),function(x) simfunc_stoch(input,rep=x))) # without dplyr
  sol <- optim(par=transto(p),fn=KR_val_func,input=input,dat=bootdat,control=list(maxit=2000))
  names(sol$par) <- genotype_names[input$pindex]
  paramfit <- transfrom(sol$par)
  paramfit
}