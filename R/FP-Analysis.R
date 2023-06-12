# FP_Analysis.R

# A script to analyze florescence polarization data for protein-nucleic acid binding and competition

###############

##### DEVELOPER NOTES

# None

###############

FPalyze <- function(experiment.type,
                    path.to.file='./',
                    file.name=c('par.txt','perp.txt'),
                    enzyme='Predator',
                    prey.molecule='Prey',
                    decoy.molecule='Decoy',
                    parameters.file=NULL,
                    save.data=T,
                    save.name='Results',
                    plot.pdf=F,
                    plot.name='Graphs',
                    save.console=T,
                    save.id='Console',
                    variant.concentrations=NULL,
                    data.size=NULL,
                    time.step=30,
                    t.zero=1,
                    incubation.time=NULL,
                    coerce.timepoints=F,
                    equilibrium.points=10,
                    outliers='none',
                    scale.data=F,
                    default.mP.values=F,
                    background.subtraction=F,
                    G.factor=1,
                    reverse.timepoints=F,
                    regress.data=T,
                    coerce.offrates=T,
                    match.optimizers=T,
                    manual.fEP.adjust=F,
                    FP.baseline=NULL,
                    legend.location='bottomright',
                    total.P=NULL,
                    show.constants=F,
                    use.anisotropy=T,
                    exponentials=1,
                    fit.DT.mdls=T,
                    drift.corr=F){

  DED.fit <- function(data){
    FXNa = FP ~ Mn + (1 - Mn)*s*exp(-lambda.1*tt) + (1 - Mn)*(1 - s)*exp(-lambda.1*beta*tt)
    FXNb = FP ~ 0.5*((exp(-1*lambda.1*tt)+Mn1*(1-exp(-1*lambda.1*tt))) + (exp(-1*lambda.1*beta*tt)+Mn2*(1-exp(-1*lambda.1*beta*tt))))
    # step 1a
    mdl.pre = smooth.spline(x = c(0,data$tt),y = c(1,data$FP),spar = 0.5)
    minAVG <- min(predict(mdl.pre,x=data$tt,deriv = 0)$y)
    # step 1b
    start.grid=expand.grid(lambda.1 = 10^c(-4:0),
                           beta = 2^c(-8:-1),
                           Mn = minAVG,
                           s = seq(0.1,0.9,0.1))
    mdl.1=nls2::nls2(FXNa,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-6,warnOnly=TRUE))
    # step 2
    start.grid=expand.grid(lambda.1 = coefficients(mdl.1)[['lambda.1']]*10^seq(-1,1,len = 11),
                           beta = coefficients(mdl.1)[['beta']]*2^seq(-1,1,len = 11),
                           Mn = coefficients(mdl.1)[['Mn']]+seq(-0.05,0.05,len = 11),
                           s = coefficients(mdl.1)[['s']]+seq(-0.05,0.05,len = 11))
    mdl.2=nls2::nls2(FXNa,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-6,warnOnly=TRUE))
    # step 3
    start.Mn1 = 1 - 2*(1 - coefficients(mdl.2)[['Mn']])*coefficients(mdl.2)[['s']]
    start.Mn2 = 2*coefficients(mdl.2)[['Mn']] + 2*(1 - coefficients(mdl.2)[['Mn']])*coefficients(mdl.2)[['s']] - 1
    start.grid=expand.grid(lambda.1 = coefficients(mdl.2)[['lambda.1']]*10^seq(-0.2,0.2,len = 6),
                           beta = coefficients(mdl.2)[['beta']]*2^seq(-0.2,0.2,len = 6),
                           Mn1 = start.Mn1+seq(-0.01,0.01,len = 6),
                           Mn2 = start.Mn2+seq(-0.01,0.01,len = 6))
    mdl.3=nls2::nls2(FXNb,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-6,warnOnly=TRUE))
    # step 4
    try(mdl.4 <- nls2::nls2(FXNb,data,start=list('lambda.1'=coefficients(mdl.3)[['lambda.1']],'beta'=coefficients(mdl.3)[['beta']],'Mn1'=coefficients(mdl.3)[['Mn1']],'Mn2'=coefficients(mdl.3)[['Mn2']]),control=list(minFactor=1e-10,maxiter=1e3,warnOnly=TRUE)),silent = TRUE)
    if (exists('mdl.4')==T){
      return(mdl.4)
    } else {
      return(mdl.3)
    }
  }

  DT.fit <- function(data,match.optimizers){
    FXNnull = y~k.n1*(x*1e-9)^N/((x*1e-9)^N+(Bhalf*1e-9)^N)
    FXNdt = y~k.n1*(x*1e-9)^N/((x*1e-9)^N+(Bhalf*1e-9)^N)+k.theta*(x*1e-9)
    # step 1a
    start.grid=expand.grid(k.n1 = 10^c(-5:0),
                           k.theta = 10^c(0:4),
                           N = 2^c(-3:3),
                           Bhalf = 1e4*4^c(-2:-5))
    mdl.DT1=nls2::nls2(FXNdt,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-10,warnOnly=TRUE))
    # step 1b
    start.grid=expand.grid(k.n1 = coefficients(mdl.DT1)[['k.n1']]*10^seq(-1,1,len=6),
                           k.theta = coefficients(mdl.DT1)[['k.theta']]*10^seq(-1,1,len=6),
                           N = coefficients(mdl.DT1)[['N']]*2^seq(-1,1,len=6),
                           Bhalf = coefficients(mdl.DT1)[['Bhalf']]*4^seq(-1,1,len=6))
    mdl.DT2=nls2::nls2(FXNdt,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-10,warnOnly=TRUE))
    # step 1c
    start.grid=expand.grid(k.n1 = coefficients(mdl.DT2)[['k.n1']]*10^seq(-0.4,0.4,len=6),
                           k.theta = coefficients(mdl.DT2)[['k.theta']]*10^seq(-0.4,0.4,len=6),
                           N = coefficients(mdl.DT2)[['N']]*2^seq(-0.4,0.4,len=6),
                           Bhalf = coefficients(mdl.DT2)[['Bhalf']]*4^seq(-0.4,0.4,len=6))
    mdl.DT3=nls2::nls2(FXNdt,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-10,warnOnly=TRUE))
    # step 1d
    try(mdl.DT4 <- nls2::nls2(FXNdt,data,start=list('k.n1'=coefficients(mdl.DT3)[['k.n1']],'k.theta'=coefficients(mdl.DT3)[['k.theta']],'N'=coefficients(mdl.DT3)[['N']],'Bhalf'=coefficients(mdl.DT3)[['Bhalf']]),control=list(minFactor=1e-10,maxiter=1e3,warnOnly=TRUE),lower=c(0,0,0.1,1),upper=c(1,1e4,10,2.5e3),algorithm = 'port'),silent = TRUE)
    if (exists('mdl.DT4')==T){
      DT.mdl = mdl.DT4
    } else {
      DT.mdl = mdl.DT3
    }
    # step 2a
    if(match.optimizers==T){
      data=list('x'=data$x,'y'=data$y,'N'=coefficients(DT.mdl)[['N']],'Bhalf'=coefficients(DT.mdl)[['Bhalf']])
      start.grid=expand.grid(k.n1 = 10^c(-5:0))
      mdl.NULL1=nls2::nls2(FXNnull,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-10,warnOnly=TRUE))
    }
    if(match.optimizers==F){
      start.grid=expand.grid(k.n1 = 10^c(-5:0),
                             N = 2^c(-3:3),
                             Bhalf = 1e4*4^c(-2:-5))
      mdl.NULL1=nls2::nls2(FXNnull,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-10,warnOnly=TRUE))
    }
    # step 2b
    if(match.optimizers==T){
      start.grid=expand.grid(k.n1 = coefficients(mdl.NULL1)[['k.n1']]*10^seq(-1,1,len=6))
      mdl.NULL2=nls2::nls2(FXNnull,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-10,warnOnly=TRUE))
    }
    if(match.optimizers==F){
      start.grid=expand.grid(k.n1 = coefficients(mdl.NULL1)[['k.n1']]*10^seq(-1,1,len=6),
                             N = coefficients(mdl.NULL1)[['N']]*2^seq(-1,1,len=6),
                             Bhalf = coefficients(mdl.NULL1)[['Bhalf']]*4^seq(-1,1,len=6))
      mdl.NULL2=nls2::nls2(FXNnull,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-10,warnOnly=TRUE))
    }
    # step 2c
    if(match.optimizers==T){
      start.grid=expand.grid(k.n1 = coefficients(mdl.NULL2)[['k.n1']]*10^seq(-0.4,0.4,len=6))
      mdl.NULL3=nls2::nls2(FXNnull,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-10,warnOnly=TRUE))
    }
    if(match.optimizers==F){
      start.grid=expand.grid(k.n1 = coefficients(mdl.NULL2)[['k.n1']]*10^seq(-0.4,0.4,len=6),
                             N = coefficients(mdl.NULL2)[['N']]*2^seq(-0.4,0.4,len=6),
                             Bhalf = coefficients(mdl.NULL2)[['Bhalf']]*4^seq(-0.4,0.4,len=6))
      mdl.NULL3=nls2::nls2(FXNnull,data,start = start.grid,algorithm = 'brute-force',control=list(minFactor=1e-10,warnOnly=TRUE))
    }
    # step 2d
    if(match.optimizers==T){
      try(mdl.NULL4 <- nls2::nls2(FXNnull,data,start=list('k.n1'=coefficients(mdl.NULL3)[['k.n1']]),control=list(minFactor=1e-10,maxiter=1e3,warnOnly=TRUE),lower=c(0),upper=c(1),algorithm = 'port'),silent = TRUE)
    }
    if(match.optimizers==F){
      try(mdl.NULL4 <- nls2::nls2(FXNnull,data,start=list('k.n1'=coefficients(mdl.NULL3)[['k.n1']],'N'=coefficients(mdl.NULL3)[['N']],'Bhalf'=coefficients(mdl.NULL3)[['Bhalf']]),control=list(minFactor=1e-10,maxiter=1e3,warnOnly=TRUE),lower=c(0,0.1,1),upper=c(1,10,2.5e3),algorithm = 'port'),silent = TRUE)
    }
    if (exists('mdl.DT4')==T){
      NULL.mdl = mdl.NULL4
    } else {
      NULL.mdl = mdl.NULL3
    }
    # step 3
    MDLs=list('NULL'=NULL.mdl,'DT'=DT.mdl)
    return(MDLs)
  }

  IC50.fit <- function(data){
    FXN <- y ~ 1 - (1-low)*(D^h)/(D^h + IC50^h)
    start.grid=expand.grid(IC50 = 1e4*2^c(0:-11),
                           h = 2^seq(-2,2,len = 10),
                           low = seq(-1,1,0.1))
    try(mdl <- nls2::nls2(FXN,data,start=data.frame('IC50'=start.grid$IC50,'h'=start.grid$h,'low'=start.grid$low),control=list(minFactor=1e-10,maxiter=1e3,warnOnly=TRUE)),silent = TRUE)
    if (exists('mdl')==T){
      return(mdl)
    } else {
      return(NULL)
    }
  }

  BIND.fit <- function(data,total.P,experiment.type){
    FXNstd <- FP~(E/(E+Kd))*(Mx-Mn)+Mn
    FXNhill <- FP~(E^n/(E^n+Kd))*(Mx-Mn)+Mn
    FXNquad <- FP~(Mx-Mn)*(((Pt+(E+Kd))-((Pt+(E+Kd))^2-4*Pt*E)^0.5)/(2*Pt))+Mn
    FXNstoich <- FP~(Mx-Mn)*(((Pt+S*(E+Kd))-((Pt+S*(E+Kd))^2-4*Pt*S*E)^0.5)/(2*Pt))+Mn
    # STD
    if(experiment.type=='Kd'){
      start.data=(list('FP'=data$FP,'E'=data$E))
      start.grid=expand.grid(Kd = 1e3*2^c(3:-11),
                             Mx = max(data$FP,na.rm = T)*2^seq(0,3,len = 5),
                             Mn = min(data$FP,na.rm = T))
      try(std.mdl <- nls2::nls2(FXNstd,start.data,start=data.frame('Kd'=start.grid$Kd,'Mx'=start.grid$Mx,'Mn'=start.grid$Mn),control=list(minFactor=1e-10,maxiter=1e3,warnOnly=TRUE)),silent = TRUE)
    }
    # HILL
    if(experiment.type=='Kd'){
      start.data=(list('FP'=data$FP,'E'=data$E))
      start.grid=expand.grid(Kd = 1e3*2^c(3:-11),
                             n = 2^seq(-2,2,len = 10),
                             Mx = max(data$FP,na.rm = T)*2^seq(0,3,len = 5),
                             Mn = min(data$FP,na.rm = T))
      try(hill.mdl <- nls2::nls2(FXNhill,start.data,start=data.frame('Kd'=start.grid$Kd,'n'=start.grid$n,'Mx'=start.grid$Mx,'Mn'=start.grid$Mn),control=list(minFactor=1e-10,maxiter=1e3,warnOnly=TRUE)),silent = TRUE)
    }
    # QUAD
    if(experiment.type=='Kd' & is.null(total.P)==F){
      start.grid=expand.grid(Kd = 1e3*2^c(3:-11),
                             Mx = max(data$FP,na.rm = T)*2^seq(0,3,len = 5),
                             Mn = min(data$FP,na.rm = T))
      try(quad.mdl <- nls2::nls2(FXNquad,data,start=data.frame('Kd'=start.grid$Kd,'Mx'=start.grid$Mx,'Mn'=start.grid$Mn),control=list(minFactor=1e-10,maxiter=1e3,warnOnly=TRUE)),silent = TRUE)
    }
    # STOICH
    if(experiment.type=='STOICH' & is.null(total.P)==F){
      start.grid=expand.grid(Kd = 1e3*2^c(3:-11),
                             S = 2^seq(-2,2,len = 10),
                             Mx = max(data$FP,na.rm = T)*2^seq(0,3,len = 5),
                             Mn = min(data$FP,na.rm = T))
      try(stoich.mdl <- nls2::nls2(FXNstoich,data,start=data.frame('Kd'=start.grid$Kd,'S'=start.grid$S,'Mx'=start.grid$Mx,'Mn'=start.grid$Mn),control=list(minFactor=1e-10,maxiter=1e3,warnOnly=TRUE)),silent = TRUE)
    }
    # Reporting
    if(exists('std.mdl')==F){
      std.mdl=NULL
      if(experiment.type=='Kd'){
        warning(paste0('Binding curve regression failed with Standard Model.'))
      }
    }
    if(exists('hill.mdl')==F){
      hill.mdl=NULL
      if(experiment.type=='Kd'){
        warning(paste0('Binding curve regression failed with Hill Model.'))
      }
    }
    if(exists('quad.mdl')==F){
      quad.mdl=NULL
      if(experiment.type=='Kd' & is.null(total.P)==F){
        warning(paste0('Binding curve regression failed with Quadratic Model.'))
      }
    }
    if(exists('stoich.mdl')==F){
      stoich.mdl=NULL
      if(experiment.type=='STOICH' & is.null(total.P)==F){
        warning(paste0('Binding curve regression failed with Stoichiometry Model.'))
      }
      if(experiment.type=='STOICH' & is.null(total.P)==T){
        warning(paste0('Regression with Stoichiometry Model not possible -- please provide ligand concentration (total.P).'))
      }
    }
    BIND.mdls=list('STD'=std.mdl,'HILL'=hill.mdl,'QUAD'=quad.mdl,'STOICH'=stoich.mdl)
    return(BIND.mdls)
  }

  ##### BEGIN SCRIPT

  if (is.null(parameters.file)==F){

    source(file = paste0(path.to.file,parameters.file))

  }

  if (is.null(parameters.file)==T){

    if ((experiment.type=='Kd' | experiment.type=='STOICH') & is.null(variant.concentrations)==T){
      variant.concentrations=signif(c(1e3/2^c(0:10),0),3) # what are the concentrations of the variant molecule, nM
    }
    if (experiment.type=='COMP' & is.null(variant.concentrations)==T){
      variant.concentrations=signif(c(1e4/2^c(0:10),0),3) # what are the concentrations of the variant molecule, nM
    }

  }

  ### Data Import and Analysis

  # Import and clean raw data
  if (default.mP.values==T){
    if (coerce.timepoints==F){
      raw=as.data.frame(read.csv(file = paste(path.to.file,file.name,sep = ''),header = TRUE,sep = '\t'))
    }
    if (coerce.timepoints==T){
      t=seq(t.zero*60,t.zero*60+incubation.time*60,time.step)
      raw=as.data.frame(read.csv(file = paste(path.to.file,file.name,sep = ''),header = TRUE,sep = '\t'))[1:length(t),]
    }
    if (outliers[1]!='none'){
        raw[,colnames(raw) %in% outliers]=NA
    }
    if(is.null(data.size)==T){
      if(ncol(raw)==13){
        data.size='single'
      }
      if(ncol(raw)==25){
        data.size='half'
      }
      if(ncol(raw)==49){
        data.size='full'
      }
    }

    # Generate time vector
    if(is.null(incubation.time)==T){
      incubation.time=time.step/60*(nrow(raw)-1)
      t=seq(t.zero*60,t.zero*60+incubation.time*60,time.step)
    }

    # Import and clean raw data, continued
    if (data.size=='full'){
      data=array(0,dim = c(length(t),12,4))
    }
    if (data.size=='half'){
      data=array(0,dim = c(length(t),12,2))
    }
    if (data.size=='single'){
      data=array(0,dim = c(length(t),12,1))
    }
    if (data.size=='single'){
      data[,,1]=as.matrix(raw[,1+seq(1,12,1)])
    }
    if (data.size=='full' | data.size=='half'){
      data[,,1]=as.matrix(raw[,1+seq(1,24,2)])
      data[,,2]=as.matrix(raw[,1+seq(2,24,2)])
    }
    if (data.size=='full'){
      data[,,3]=as.matrix(raw[,1+seq(25,48,2)])
      data[,,4]=as.matrix(raw[,1+seq(26,48,2)])
    }
    if(length(table(data=='Invalid'))>1){
      data[data=='Invalid']=NA; data=array(as.numeric(data),dim = dim(data))
    }
  }
  if (default.mP.values==F){
    if (coerce.timepoints==F){
      raw.par=as.data.frame(read.csv(file = paste(path.to.file,file.name[1],sep = ''),header = TRUE,sep = '\t'))
      raw.perp=as.data.frame(read.csv(file = paste(path.to.file,file.name[2],sep = ''),header = TRUE,sep = '\t'))
    }
    if (coerce.timepoints==T){
      if(is.null(incubation.time)==T){
        stop('Coercion of time points not possible -- please manually specificy time parameters (incubation.time, time.step)')
      }
      t=seq(t.zero*60,t.zero*60+incubation.time*60,time.step)
      raw.par=as.data.frame(read.csv(file = paste(path.to.file,file.name[1],sep = ''),header = TRUE,sep = '\t'))[1:length(t),]
      raw.perp=as.data.frame(read.csv(file = paste(path.to.file,file.name[2],sep = ''),header = TRUE,sep = '\t'))[1:length(t),]
    }
    if (outliers[1]!='none'){
        raw.par[,colnames(raw.par) %in% outliers]=NA
        raw.perp[,colnames(raw.perp) %in% outliers]=NA
    }
    if(is.null(data.size)==T){
      if(ncol(raw.par)==13){
        data.size='single'
      }
      if(ncol(raw.par)==25){
        data.size='half'
      }
      if(ncol(raw.par)==49){
        data.size='full'
      }
    }

    # Generate time vector
    if(is.null(incubation.time)==T){
      incubation.time=time.step/60*(nrow(raw.par)-1)
      t=seq(t.zero*60,t.zero*60+incubation.time*60,time.step)
    }

    # Import and clean raw data, continued
    if (data.size==('full')){
      data.par=array(0,dim = c(length(t),12,4))
    }
    if (data.size==('half')){
      data.par=array(0,dim = c(length(t),12,2))
    }
    if (data.size==('single')){
      data.par=array(0,dim = c(length(t),12,1))
      data.par[,,1]=as.matrix(raw.par[,2:13])
    }
    if (data.size=='full' | data.size=='half'){
      data.par[,,1]=as.matrix(raw.par[,1+seq(1,24,2)])
      data.par[,,2]=as.matrix(raw.par[,1+seq(2,24,2)])
    }
    if (data.size==('full')){
      data.par[,,3]=as.matrix(raw.par[,1+seq(25,48,2)])
      data.par[,,4]=as.matrix(raw.par[,1+seq(26,48,2)])
    }
    if (data.size==('full')){
      data.perp=array(0,dim = c(length(t),12,4))
    }
    if (data.size==('half')){
      data.perp=array(0,dim = c(length(t),12,2))
    }
    if (data.size==('single')){
      data.perp=array(0,dim = c(length(t),12,1))
      data.perp[,,1]=as.matrix(raw.perp[,2:13])
    }
    if (data.size=='full' | data.size=='half'){
      data.perp[,,1]=as.matrix(raw.perp[,1+seq(1,24,2)])
      data.perp[,,2]=as.matrix(raw.perp[,1+seq(2,24,2)])
    }
    if (data.size==('full')){
      data.perp[,,3]=as.matrix(raw.perp[,1+seq(25,48,2)])
      data.perp[,,4]=as.matrix(raw.perp[,1+seq(26,48,2)])
    }
    if (((length(table((data.perp=='OVER')))>1) | (length(table((data.par=='OVER')))>1))){
      data.perp[data.perp=='OVER']=NA; data.perp=array(as.numeric(data.perp),dim = dim(data.perp))
      data.par[data.par=='OVER']=NA; data.par=array(as.numeric(data.par),dim = dim(data.par))
    }
    if (background.subtraction==T & data.size==('full')){
      data.par=data.par-mean(as.matrix(raw.par[,50:53]))
      data.perp=data.perp-mean(as.matrix(raw.perp[,50:53]))
    }
    if (background.subtraction==T & data.size==('half')){
      data.par=data.par-mean(as.matrix(raw.par[,26:29]))
      data.perp=data.perp-mean(as.matrix(raw.perp[,26:29]))
    }
    if (background.subtraction==T & data.size==('single')){
      data.par=data.par-mean(as.matrix(raw.par[,14:17]))
      data.perp=data.perp-mean(as.matrix(raw.perp[,14:17]))
    }
    data.mP=1000*(data.par-G.factor*data.perp)/(data.par+G.factor*data.perp)
    data.A=(data.par-G.factor*data.perp)/(data.par+2*G.factor*data.perp)
    if(use.anisotropy==F){
      data=data.mP
    }
    if(use.anisotropy==T){
      data=data.A
    }
  }

  # Data drift correction
  if (drift.corr==T & experiment.type=='COMP'){
    ref.a=rowMeans(data[,variant.concentrations==0,])
    ref.b=ref.a[length(ref.a)]-ref.a
    ref.c=array(ref.b,dim = dim(data))
    data=data+ref.c
    rm(ref.a,ref.b,ref.c)
  }

  # Scale raw data
  if (experiment.type=='Kd' | experiment.type=='STOICH'){
    if (scale.data==T){
      if (reverse.timepoints==F){
        data.scaled=(data-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.min(variant.concentrations),]))))/(mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.min(variant.concentrations),]))))
      }
      if (reverse.timepoints==T){
        data.scaled=(data-mean(na.omit(c(data[1:3,which.min(variant.concentrations),]))))/(mean(na.omit(c(data[1:3,which.max(variant.concentrations),])))-mean(na.omit(c(data[1:3,which.min(variant.concentrations),]))))
      }
    }
    if (scale.data==F){
        data.scaled=data
    }
  }
  if (experiment.type=='COMP'){
    data.scaled=(data-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))/(mean(na.omit(c(data[1,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))
  }

  if (experiment.type=='Kd' | experiment.type=='STOICH'){
    # Calculate equilibrium values from scaled data
    if (reverse.timepoints==F & data.size!='single'){
      data.scaled.eq=apply(data.scaled[(length(t)-equilibrium.points+1):length(t),,],c(2,3),mean)
    }
    if (reverse.timepoints==T & data.size!='single'){
      data.scaled.eq=apply(data.scaled[1:equilibrium.points,,],c(2,3),mean)
    }
    if (reverse.timepoints==F & data.size=='single'){
      data.scaled.eq=apply(data.scaled[(length(t)-equilibrium.points+1):length(t),,],c(2),mean)
    }
    if (reverse.timepoints==T & data.size=='single'){
      data.scaled.eq=apply(data.scaled[1:equilibrium.points,,],c(2),mean)
    }

    # Fit equilibrium binding curve data
    if (experiment.type=='STOICH' | experiment.type=='Kd'){
      if (data.size=='full' | data.size=='half'){model.data=list('E'=rep(variant.concentrations,times=dim(data)[3])[is.na(c(data.scaled.eq))==FALSE],'FP'=c(data.scaled.eq)[is.na(c(data.scaled.eq))==FALSE],'Pt'=total.P)}
      if (data.size=='single'){model.data=list('E'=rep(variant.concentrations,times=1)[is.na(c(data.scaled.eq))==FALSE],'FP'=c(data.scaled.eq)[is.na(c(data.scaled.eq))==FALSE],'Pt'=total.P)}
    }

    suppressWarnings(BIND.models <- BIND.fit(model.data,total.P,experiment.type),classes = 'all')

  }

  ### Report Results

  # Open txt file to save console output
  if (save.console==T){
    sink(paste0(path.to.file,save.id,'.txt'))
  }

  # Open pdf to save plots
  if (plot.pdf==TRUE){
    pdf(file = paste(path.to.file,plot.name,'.pdf',sep = ''))
  }

  if (experiment.type=='Kd'){
    # Plot association curves
    par(mfrow=c(4,3))
    for (i in 1:12){
      if (scale.data==T){
        plot(rep(t/60,if(data.size!='single'){times=dim(data)[3]}else{times=1}),c(data.scaled[,i,]),type = 'p',main = paste(enzyme,'-',prey.molecule,' Association Curve:  ',enzyme,' at ',variant.concentrations[i],' nM',sep = ''),xlab = 'Time (minutes)',ylab = if(use.anisotropy==T){'Relative Anisotropy'}else{'Relative Polarization'},ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex.main=0.8,cex=1,col='black')
      }
      if (scale.data==F){
        plot(rep(t/60,if(data.size!='single'){times=dim(data)[3]}else{times=1}),c(data.scaled[,i,]),type = 'p',main = paste(enzyme,'-',prey.molecule,' Association Curve:  ',enzyme,' at ',variant.concentrations[i],' nM',sep = ''),xlab = 'Time (minutes)',ylab = if(use.anisotropy==T){'Anisotropy'}else{'Polarization (mP)'},ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex.main=0.8,cex=1,col='black')
      }
      lines(t/60,rowMeans(data.scaled[,i,],na.rm = TRUE),col='red',lwd=4)
    }

    # Plot equilibrium binding curve
    par(mfrow=c(2,1))
    if (scale.data==T){
      plot(rep(variant.concentrations,if(data.size!='single'){times=dim(data)[3]}else{times=1}),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = if(use.anisotropy==T){'Relative Anisotropy'}else{'Relative Polarization'},ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
    }
    if (scale.data==F){
      plot(rep(variant.concentrations,if(data.size!='single'){times=dim(data)[3]}else{times=1}),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = if(use.anisotropy==T){'Anisotropy'}else{'Polarization (mP)'},ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
    }
    try(lines(min(variant.concentrations):max(variant.concentrations),predict(BIND.models[['STD']],newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lwd=3),silent = TRUE)
    try(lines(min(variant.concentrations):max(variant.concentrations),predict(BIND.models[['HILL']],newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='red',lwd=3),silent = TRUE)
    try(lines(min(variant.concentrations):max(variant.concentrations),predict(BIND.models[['QUAD']],newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='green',lwd=3),silent = TRUE)
    legend('bottomright',legend = c('Std.','Hill','Quad.'),fill = c('blue','red','green'),col = c('blue','red','green'),cex=1.5)

    if (scale.data==T){
      plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = if(use.anisotropy==T){'Relative Anisotropy'}else{'Relative Polarization'},ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
    }
    if (scale.data==F){
      plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = if(use.anisotropy==T){'Anisotropy'}else{'Polarization (mP)'},ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
    }
    arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
    try(lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(BIND.models[['STD']],newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lwd=3),silent = TRUE)
    try(lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(BIND.models[['HILL']],newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='red',lwd=3),silent = TRUE)
    try(lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(BIND.models[['QUAD']],newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='green',lwd=3),silent = TRUE)
    legend('topleft',legend = c('Std.','Hill','Quad.'),fill = c('blue','red','green'),col = c('blue','red','green'),cex=1.5)

    # Summarize binding curve regression
    print(summary(BIND.models[['STD']])); print(paste(rep('#',times=100),collapse = ''),quote = F)
    print(summary(BIND.models[['HILL']])); print(paste(rep('#',times=100),collapse = ''),quote = F)
    print(summary(BIND.models[['QUAD']])); print(paste(rep('#',times=100),collapse = ''),quote = F)

    # Model comparison statistics
    if(is.null(BIND.models[['STD']])==F & is.null(BIND.models[['HILL']])==F & (is.null(BIND.models[['QUAD']])==F | is.null(total.P)==T)){
      BIC=BIC(BIND.models[['STD']],BIND.models[['HILL']])
      print(BIC,quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      delta.BIC=BIC(BIND.models[['HILL']])-BIC(BIND.models[['STD']])
      print(paste('ΔBIC = ',delta.BIC,sep = ''),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      if (delta.BIC<=-10){
        print('These data favor the HILL model!',quote = F)
        print(paste0('Kd = ',signif(summary(BIND.models[['HILL']])[['coefficients']][1,1],2),' ± ',signif(summary(BIND.models[['HILL']])[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste0('n = ',signif(summary(BIND.models[['HILL']])[['coefficients']][2,1],2),' ± ',signif(summary(BIND.models[['HILL']])[['coefficients']][2,2],2)),quote=F)
        print(paste0('B_0.5 ≈ ',signif(summary(BIND.models[['HILL']])[['coefficients']][1,1]^(1/summary(BIND.models[['HILL']])[['coefficients']][2,1]),2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (delta.BIC>=5){
        print('These data favor the STANDARD model!',quote = F)
        print(paste0('Kd = ',signif(summary(BIND.models[['STD']])[['coefficients']][1,1],2),' ± ',signif(summary(BIND.models[['STD']])[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (delta.BIC>-10 & delta.BIC<5){
        print('EITHER model may have produced these data!',quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        print(paste0('Kd = ',signif(summary(BIND.models[['STD']])[['coefficients']][1,1],2),' ± ',signif(summary(BIND.models[['STD']])[['coefficients']][1,2],2),' nM'),quote=F)
        print('OR',quote=F)
        print(paste0('Kd = ',signif(summary(BIND.models[['HILL']])[['coefficients']][1,1],2),' ± ',signif(summary(BIND.models[['HILL']])[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste0('n = ',signif(summary(BIND.models[['HILL']])[['coefficients']][2,1],2),' ± ',signif(summary(BIND.models[['HILL']])[['coefficients']][2,2],2)),quote=F)
        print(paste0('B_0.5 ≈ ',signif(summary(BIND.models[['HILL']])[['coefficients']][1,1]^(1/summary(BIND.models[['HILL']])[['coefficients']][2,1]),2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
    }

  }

  if (experiment.type=='COMP'){
    if (exponentials==1){
    # Plot dissociation curves
    par(mfrow=c(4,3),mar=c(4.5,5,3,1),oma = c(0,0,5,0))
    cmp.models=list(NULL)
    cmp.models.coefficients.lambda=matrix(NA,ncol=dim(data)[3],nrow=length(variant.concentrations)); cmp.models.coefficients.Mn=matrix(NA,if(data.size!='single'){ncol=dim(data)[3]}else{ncol=1},nrow=length(variant.concentrations))
    for (i in 1:12){
      if (scale.data==T){
        plot(rep(t/60,times=dim(data)[3]),c(data.scaled[,i,]),type = 'p',main = paste0('[Decoy] = ',signif(variant.concentrations[i]/1e3,2),' µM'),xlab = 'Time (min)',ylab = if(use.anisotropy==T){'Relative Anisotropy'}else{'Relative Polarization'},ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex=1,cex.main=2.5,cex.axis=2,cex.lab=2)
      }
      if (scale.data==F){
        plot(rep(t/60,times=dim(data)[3]),c(data[,i,]),type = 'p',main = paste0('[Decoy] = ',signif(variant.concentrations[i]/1e3,2),' µM'),xlab = 'Time (min)',ylab = if(use.anisotropy==T){'Anisotropy'}else{'Polarization (mP)'},ylim = c(min(na.omit(c(data))),max(na.omit(c(data)))),cex=1,cex.main=2.5,cex.axis=2,cex.lab=2)
      }
      if (i==1){
        title(paste0(enzyme,':  ',prey.molecule,' --> ',decoy.molecule),line = 1.5,outer = TRUE,cex.main = 4)
      }
      if (regress.data==T){
        for (q in 1:dim(data)[3]){
          cmp.models.data=list('FP'=c(data.scaled[,i,q]),'tt'=t); try(cmp.models[[paste(i,q,sep = '_')]]<-nls(FP~(1-Mn)*exp(-kn1*tt)+Mn,data=cmp.models.data,start=list('kn1'=1e-3,'Mn'=0),control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE)),silent = TRUE)
          if (is.null(cmp.models[[paste(i,q,sep = '_')]])==F){
            cmp.models.coefficients.lambda[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[1]; cmp.models.coefficients.Mn[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[2]
          }
          if (is.null(cmp.models[[paste(i,q,sep = '_')]])==T){
            cmp.models.coefficients.lambda[i,q]=0; cmp.models.coefficients.Mn[i,q]=1
            warning(paste0('Reaction ',i,'-',q,' does not fit an exponential decay -- it was coerced to a horizontal line at binding saturation.'))
          }
        }
        cmp.models.data=list('FP'=c(data.scaled[,i,]),'tt'=rep(t,times=dim(data)[3])); try(cmp.models[[paste(i,'all',sep = '_')]]<-nls(FP~(1-Mn)*exp(-kn1*tt)+Mn,data=cmp.models.data,start=list('kn1'=1e-3,'Mn'=0),control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE)),silent = TRUE)
        if (scale.data==T & is.null(cmp.models[[paste(i,'all',sep = '_')]])==F){
          lines(t/60,predict(cmp.models[[paste(i,'all',sep = '_')]],newdata=list('tt'=t)),col='red',lwd=3)
        }
        if (scale.data==F & is.null(cmp.models[[paste(i,'all',sep = '_')]])==F){
          lines(t/60,predict(cmp.models[[paste(i,'all',sep = '_')]],newdata=list('tt'=t))*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))),col='red',lwd=3)
        }
        if (scale.data==T & is.null(cmp.models[[paste(i,'all',sep = '_')]])==T){
          lines(t/60,rep(1,times=length(t)),col='red',lwd=3)
        }
        if (scale.data==F & is.null(cmp.models[[paste(i,'all',sep = '_')]])==T){
          lines(t/60,rep(1,times=length(t))*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))),col='red',lwd=3)
        }
      }
    }

    # Plot comparative dissociation curves
    if (regress.data==T){
      par(fig=c(0,0.5,0.5,0.9),mar=c(4.5,5,3,1),oma = c(0,0,0,0))
      if (scale.data==T){
        plot(NULL,NULL,main = 'Comparison of Competition Reactions',xlab = 'Time (min)',ylab = if(use.anisotropy==T){'Relative Anisotropy'}else{'Relative Polarization'},ylim = c(0,1),xlim = c(min(t)/60,max(t)/60),cex.main=2,cex.axis=2,cex.lab=2)
      }
      if (scale.data==F){
        plot(NULL,NULL,main = 'Comparison of Competition Reactions',xlab = 'Time (min)',ylab = if(use.anisotropy==T){'Anisotropy'}else{'Polarization (mP)'},ylim = c(min(na.omit(c(data))),max(na.omit(c(data)))),xlim = c(min(t)/60,max(t)/60),cex.main=2,cex.axis=2,cex.lab=2)
      }
      for (i in 12:1){
        if (scale.data==T & is.null(cmp.models[[paste(i,'all',sep = '_')]])==F){
          lines(t/60,predict(cmp.models[[paste(i,'all',sep = '_')]],newdata=list('tt'=t)),lwd=3,col=grey.colors(12)[i])
        }
        if (scale.data==F & is.null(cmp.models[[paste(i,'all',sep = '_')]])==F){
          lines(t/60,predict(cmp.models[[paste(i,'all',sep = '_')]],newdata=list('tt'=t))*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))),lwd=3,col=grey.colors(12)[i])
        }
        if (scale.data==T & is.null(cmp.models[[paste(i,'all',sep = '_')]])==T){
          lines(t/60,rep(1,times=length(t)),lwd=3,col=grey.colors(12)[i])
        }
        if (scale.data==F & is.null(cmp.models[[paste(i,'all',sep = '_')]])==T){
          lines(t/60,rep(1,times=length(t))*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))),lwd=3,col=grey.colors(12)[i])
        }
      }
      if (scale.data==T){
        fields::colorbar.plot(x = 0.75*incubation.time,y = 0.9,strip = seq(0.3,0.9,0.6/11),col = grey.colors(12),strip.width = 0.05,strip.length = 0.35)
        text((0.75-0.35/2)*incubation.time,1,labels = paste0('[Decoy] (µM):'),adj = c(0.25,0.5),cex=2)
        text((0.75-0.35/2)*incubation.time,0.9,labels = paste0(round(max(variant.concentrations)/1e3)),adj = c(1.5,0.5),cex=2)
        text((0.75+0.35/2)*incubation.time,0.9,labels = paste0(round(min(variant.concentrations)/1e3)),adj = c(-1,0.5),cex=2)
      }
      if (scale.data==F){
        fields::colorbar.plot(x = 0.75*incubation.time,y = 0.9*diff(range(na.omit(c(data))))+min(na.omit(c(data))),strip = seq(0.3,0.9,0.6/11),col = grey.colors(12),strip.width = 0.05,strip.length = 0.35)
        text((0.75-0.35/2)*incubation.time,1*diff(range(na.omit(c(data))))+min(na.omit(c(data))),labels = paste0('[Decoy] (µM):'),adj = c(0.25,0.5),cex=2)
        text((0.75-0.35/2)*incubation.time,0.9*diff(range(na.omit(c(data))))+min(na.omit(c(data))),labels = paste0(round(max(variant.concentrations)/1e3)),adj = c(1.5,0.5),cex=2)
        text((0.75+0.35/2)*incubation.time,0.9*diff(range(na.omit(c(data))))+min(na.omit(c(data))),labels = paste0(round(min(variant.concentrations)/1e3)),adj = c(-1,0.5),cex=2)
      }
      title(paste0(enzyme,':  ',prey.molecule,' --> ',decoy.molecule),line = -5,outer = TRUE,cex.main = 4.5)

      # Plot decoy-dependence curves
      if (manual.fEP.adjust==T){
        cmp.models.coefficients.Mn=(cmp.models.coefficients.Mn*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),])))-FP.baseline)/(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-FP.baseline)
      }
      k.off=cmp.models.coefficients.lambda*(1-cmp.models.coefficients.Mn); if (coerce.offrates==T){k.off[k.off<=0]=0}
      kt=min(na.omit(k.off))
      par(fig=c(0.5,1,0.7,0.9),mar=c(4.5,5,3,1),oma = c(0,0,0,0),new=TRUE)
      plot(variant.concentrations/1e3,rowMeans(cmp.models.coefficients.lambda,na.rm = TRUE),main = "Decoy-Dependence of λ",xlab = paste0('[Decoy] (µM)'),ylab = 'λ',type='p',cex.main=2,cex.lab=2,cex.axis=2)
      arrows(x0 = variant.concentrations/1e3, x1 = variant.concentrations/1e3, y0 = rowMeans(cmp.models.coefficients.lambda,na.rm = TRUE) - apply(cmp.models.coefficients.lambda,MARGIN = c(1),FUN = sd,na.rm = TRUE), y1 = rowMeans(cmp.models.coefficients.lambda,na.rm = TRUE) + apply(cmp.models.coefficients.lambda,MARGIN = c(1),FUN = sd,na.rm = TRUE),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
      par(fig=c(0.5,1,0.5,0.7),mar=c(4.5,5,3,1),oma = c(0,0,0,0),new=TRUE)
      plot(variant.concentrations/1e3,rowMeans(cmp.models.coefficients.Mn,na.rm = TRUE),main = 'Decoy-Dependence of [EP] Equilibrium',xlab = paste0('[Decoy] (µM)'),ylab = expression('Fraction [EP]'[0]),type='p',cex.main=2,cex.lab=2,cex.axis=2)
      arrows(x0 = variant.concentrations/1e3, x1 = variant.concentrations/1e3, y0 = rowMeans(cmp.models.coefficients.Mn,na.rm = TRUE) - apply(cmp.models.coefficients.Mn,MARGIN = c(1),FUN = sd,na.rm = TRUE), y1 = rowMeans(cmp.models.coefficients.Mn,na.rm = TRUE) + apply(cmp.models.coefficients.Mn,MARGIN = c(1),FUN = sd,na.rm = TRUE),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
      if (fit.DT.mdls==T){
        IC50.data=list('y'=c(cmp.models.coefficients.Mn),'D'=rep(variant.concentrations,times=ncol(cmp.models.coefficients.Mn)))
        IC50.mdl=IC50.fit(IC50.data)
        try(lines(seq(0,max(variant.concentrations)*1e-3,1e-3),predict(IC50.mdl,newdata=list('D' = seq(0,max(variant.concentrations),1))),lwd=3))
        if (show.constants==T & is.null(IC50.mdl)==F){
          cc.temp = signif(coefficients(IC50.mdl)[['IC50']],2)
          text(max(variant.concentrations)*1e-3*0.6,0.85,labels = substitute(paste('IC'[50],' = ',cc.temp,' nM',sep = "")),adj = c(0,0.5),cex=2)
          cc.temp = signif(coefficients(IC50.mdl)[['h']],2)
          text(max(variant.concentrations)*1e-3*0.6,0.6,labels = substitute(paste('h = ',cc.temp,sep = "")),adj = c(0,0.5),cex=2)
        }
      }

      par(fig=c(0,1,0,0.5),mar=c(5,6,3,1),oma = c(0,0,0,0),new=TRUE)
      plot(variant.concentrations*1e-3,rowMeans(k.off,na.rm = TRUE)*1e3,main = expression('Decoy-Dependence of k'[off]*''^obs),xlab = paste0('[Decoy] (µM)'),ylab = expression('k'['off']*''^obs*' x10'^-3*' (s'^-1*')'),type='p',cex.main=3,cex.lab=2.5,cex.axis=2.5)
      arrows(x0 = variant.concentrations*1e-3, x1 = variant.concentrations*1e-3, y0 = rowMeans(k.off,na.rm = TRUE)*1e3 - apply(k.off,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3, y1 = rowMeans(k.off,na.rm = TRUE)*1e3 + apply(k.off,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3,code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
      if (fit.DT.mdls==T){
        # Generate competition models
        DTmod.data=list('x'=(rep(variant.concentrations,times=dim(data)[3]))[is.na(c(k.off-kt))==F],'y'=na.omit(c(k.off-kt)))
        DTmod.both=DT.fit(DTmod.data,match.optimizers)
        DTmod.opt=DTmod.both[['NULL']]
        DTmod.opt2=DTmod.both[['DT']]
        lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod.opt,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt*1e3,col='green',lwd=3)
        lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod.opt2,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt*1e3,col='purple',lwd=3)
        legend(legend.location,legend = c('Classic Competition Model','Direct Transfer Model'),col = c('green','purple'),fill = c('green','purple'),cex = 2)
        if (show.constants==T){
          cc.temp = signif((summary(DTmod.opt)[['coefficients']][1,1])*1e3,2); text(max(variant.concentrations)*1e-3*0.7,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.5+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=3,col='green')
          cc.temp = signif((summary(DTmod.opt2)[['coefficients']][1,1])*1e3,2); text(max(variant.concentrations)*1e-3*0.7,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.4+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=3,col='purple')
          cc.temp = signif((summary(DTmod.opt2)[['coefficients']][4,1]),2); text(max(variant.concentrations)*1e-3*0.7,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.3+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[theta],' = ',cc.temp,' M'^-1,'s'^-1,sep = "")),adj = c(0,0.5),cex=3,col='purple')
        }

        # Report and compare models
        if (is.null(IC50.mdl)==F){
          print(paste('Equilibrium Competition Summary:',sep = ''),quote = F); print(summary(IC50.mdl),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (exists('DTmod.opt')==T){
          print(paste('Classic Model Summary:',sep = ''),quote = F); print(summary(DTmod.opt),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (exists('DTmod.opt2')==T){
          print(paste('Direct Transfer Model Summary:',sep = ''),quote = F); print(summary(DTmod.opt2),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (exists('DTmod.opt2')==T & exists('DTmod.opt')==T){
          model.comparison=BIC(DTmod.opt,DTmod.opt2)
          print(paste("Models' BIC Statistics:",sep = ''),quote = F); print(model.comparison,quote = F)
          delta.BIC=BIC(DTmod.opt2)-BIC(DTmod.opt)
          print(paste('ΔBIC = ',delta.BIC,sep = ''),quote = F)
          print(paste(rep('#',times=100),collapse = ''),quote = F)
          if (delta.BIC<=-10){
            print('These data favor the DIRECT TRANSFER model!',quote = F)
            print(paste0('k.n1 = ',signif(summary(DTmod.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print(paste0('k.theta = ',signif(summary(DTmod.opt2)[['coefficients']][4,1],2),' ± ',signif(summary(DTmod.opt2)[['coefficients']][4,2],2),' 1/M/s'),quote=F)
            print(paste(rep('#',times=100),collapse = ''),quote = F)
            if (show.constants==T){
              text(max(variant.concentrations)*1e-3*0.65,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.05+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels='BIC Favors >',cex = 5,col = 'purple',adj=c(1.5,0.5))
            }
          }
          if (delta.BIC>=5){
            print('These data favor the CLASSIC model!',quote = F)
            print(paste0('k.n1 = ',signif(summary(DTmod.opt)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print(paste(rep('#',times=100),collapse = ''),quote = F)
            if (show.constants==T){
              text(max(variant.concentrations)*1e-3*0.65,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.15+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels='BIC Favors >',cex = 5,col = 'green',adj=c(1.5,0.5))
            }
          }
          if (delta.BIC>-10 & delta.BIC<5){
            print('EITHER model may have produced these data!',quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
            print(paste0('k.n1 = ',signif(summary(DTmod.opt)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print('OR',quote=F)
            print(paste0('k.n1 = ',signif(summary(DTmod.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print(paste0('k.theta = ',signif(summary(DTmod.opt2)[['coefficients']][4,1],2),' ± ',signif(summary(DTmod.opt2)[['coefficients']][4,2],2),' 1/M/s'),quote=F)
            print(paste(rep('#',times=100),collapse = ''),quote = F)
            if (show.constants==T){
              text(max(variant.concentrations)*1e-3*0.5,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.1+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels='BIC',cex = 5,col = 'grey',adj=c(1.5,0.5))
            }
          }
        }

        # Outlier Identification
        if (outliers[1]=='none' & default.mP.values==F){
          temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(DTmod.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
          if (length(outlyers)>=1){
            print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          if (length(outlyers)==0){
            print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
        }
        if (outliers[1]=='none' & default.mP.values==T){
          temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(DTmod.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
          if (length(outlyers)>=1){
            print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          if (length(outlyers)==0){
            print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
        }
        if (outliers[1]!='none' & default.mP.values==F){
          temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(DTmod.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
          if (length(outlyers)>=1){
            print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          if (length(outlyers)==0){
            print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
        }
        if (outliers[1]!='none' & default.mP.values==T){
          temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(DTmod.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[-((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
          if (length(outlyers)>=1){
            print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          if (length(outlyers)==0){
            print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
        }
      }
    }
    }
    if (exponentials==2){
      # Plot dissociation curves
      par(mfrow=c(4,3),mar=c(4.5,5,3,1),oma = c(0,0,5,0))
      cmp.models=list(NULL)
      cmp.models.coefficients.lambda.1=matrix(NA,ncol=dim(data)[3],nrow=length(variant.concentrations))
      cmp.models.coefficients.lambda.2=cmp.models.coefficients.lambda.1
      cmp.models.coefficients.Mn.1=cmp.models.coefficients.lambda.1
      cmp.models.coefficients.Mn.2=cmp.models.coefficients.lambda.1
      for (i in 1:12){
        if (scale.data==T){
          plot(rep(t/60,times=dim(data)[3]),c(data.scaled[,i,]),type = 'p',main = paste0('[Decoy] = ',signif(variant.concentrations[i]/1e3,2),' µM'),xlab = 'Time (min)',ylab = if(use.anisotropy==T){'Relative Anisotropy'}else{'Relative Polarization'},ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex=1,cex.main=2.5,cex.axis=2,cex.lab=2)
        }
        if (scale.data==F){
          plot(rep(t/60,times=dim(data)[3]),c(data[,i,]),type = 'p',main = paste0('[Decoy] = ',signif(variant.concentrations[i]/1e3,2),' µM'),xlab = 'Time (min)',ylab = if(use.anisotropy==T){'Anisotropy'}else{'Polarization (mP)'},ylim = c(min(na.omit(c(data))),max(na.omit(c(data)))),cex=1,cex.main=2.5,cex.axis=2,cex.lab=2)
        }
        if (i==1){
          title(paste0(enzyme,':  ',prey.molecule,' --> ',decoy.molecule),line = 1.5,outer = TRUE,cex.main = 4)
        }
        if (regress.data==T){
          for (q in 1:4){
            cmp.models.data=list('FP'=c(data.scaled[,i,q]),'tt'=t)
            try(cmp.models[[paste(i,q,sep = '_')]]<-DED.fit(cmp.models.data),silent = TRUE)
            if (is.null(cmp.models[[paste(i,q,sep = '_')]])==T){
              try(cmp.models[[paste(i,q,sep = '_')]]<-nls(FP~(1-Mn)*exp(-kn1*tt)+Mn,data=cmp.models.data,start=list('kn1'=1e-3,'Mn'=0),control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE)),silent = TRUE)
              if (is.null(cmp.models[[paste(i,q,sep = '_')]])==F){
                cmp.models.coefficients.lambda.1[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[1]; cmp.models.coefficients.Mn.1[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[2]; cmp.models.coefficients.lambda.2[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[1]; cmp.models.coefficients.Mn.2[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[2]
                warning(paste0('Reaction ',i,'-',q,' does not fit a double exponential decay -- it was coerced to a single exponential.'))
              }
              if (is.null(cmp.models[[paste(i,q,sep = '_')]])==T){
                cmp.models.coefficients.lambda.1[i,q]=0; cmp.models.coefficients.Mn.1[i,q]=1; cmp.models.coefficients.lambda.1[i,q]=0; cmp.models.coefficients.Mn.2[i,q]=1; cmp.models.coefficients.lambda.2[i,q]=0
                warning(paste0('Reaction ',i,'-',q,' does not fit an exponential decay -- it was coerced to a horizontal line at binding saturation.'))
              }
            } else {
              if (is.null(cmp.models[[paste(i,q,sep = '_')]])==F){
                cmp.models.coefficients.lambda.1[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[['lambda.1']]
                cmp.models.coefficients.lambda.2[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[['lambda.1']]*coefficients(cmp.models[[paste(i,q,sep = '_')]])[['beta']]
                cmp.models.coefficients.Mn.1[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[['Mn1']]
                cmp.models.coefficients.Mn.2[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[['Mn2']]
              }
            }
          }
          cmp.models.data=list('FP'=c(data.scaled[,i,]),'tt'=rep(t,times=dim(data)[3]))
          try(cmp.models[[paste(i,'all',sep = '_')]]<-DED.fit(cmp.models.data),silent = TRUE)
          if (is.null(cmp.models[[paste(i,'all',sep = '_')]])==T){
            try(cmp.models[[paste(i,'all',sep = '_')]]<-nls(FP~(1-Mn)*exp(-kn1*tt)+Mn,data=cmp.models.data,start=list('kn1'=1e-3,'Mn'=0),control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE)),silent = TRUE)
          }
          if (scale.data==T & is.null(cmp.models[[paste(i,'all',sep = '_')]])==F){
            lines(t/60,predict(cmp.models[[paste(i,'all',sep = '_')]],newdata=list('tt'=t)),col='red',lwd=3)
          }
          if (scale.data==F & is.null(cmp.models[[paste(i,'all',sep = '_')]])==F){
            lines(t/60,predict(cmp.models[[paste(i,'all',sep = '_')]],newdata=list('tt'=t))*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))),col='red',lwd=3)
          }
        }
      }

      # Plot comparative dissociation curves
      if (regress.data==T){
        par(fig=c(0,0.5,0.5,0.9),mar=c(4.5,5,3,1),oma = c(0,0,0,0))
        if (scale.data==T){
          plot(NULL,NULL,main = 'Comparison of Competition Reactions',xlab = 'Time (min)',ylab = if(use.anisotropy==T){'Relative Anisotropy'}else{'Relative Polarization'},ylim = c(0,1),xlim = c(min(t)/60,max(t)/60),cex.main=2,cex.axis=2,cex.lab=2)
        }
        if (scale.data==F){
          plot(NULL,NULL,main = 'Comparison of Competition Reactions',xlab = 'Time (min)',ylab = if(use.anisotropy==T){'Anisotropy'}else{'Polarization (mP)'},ylim = c(min(na.omit(c(data))),max(na.omit(c(data)))),xlim = c(min(t)/60,max(t)/60),cex.main=2,cex.axis=2,cex.lab=2)
        }
        for (i in 12:1){
          if (scale.data==T & is.null(cmp.models[[paste(i,'all',sep = '_')]])==F){
            lines(t/60,predict(cmp.models[[paste(i,'all',sep = '_')]],newdata=list('tt'=t)),lwd=3,col=grey.colors(12)[i])
          }
          if (scale.data==F & is.null(cmp.models[[paste(i,'all',sep = '_')]])==F){
            lines(t/60,predict(cmp.models[[paste(i,'all',sep = '_')]],newdata=list('tt'=t))*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))),lwd=3,col=grey.colors(12)[i])
          }
          if (scale.data==T & is.null(cmp.models[[paste(i,'all',sep = '_')]])==T){
            lines(t/60,rep(1,times=length(t)),lwd=3,col=grey.colors(12)[i])
          }
          if (scale.data==F & is.null(cmp.models[[paste(i,'all',sep = '_')]])==T){
            lines(t/60,rep(1,times=length(t))*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))),lwd=3,col=grey.colors(12)[i])
          }
        }
        if (scale.data==T){
          fields::colorbar.plot(x = 0.75*incubation.time,y = 0.9,strip = seq(0.3,0.9,0.6/11),col = grey.colors(12),strip.width = 0.05,strip.length = 0.35)
          text((0.75-0.35/2)*incubation.time,1,labels = paste0('[Decoy] (µM):'),adj = c(0.25,0.5),cex=2)
          text((0.75-0.35/2)*incubation.time,0.9,labels = paste0(round(max(variant.concentrations)/1e3)),adj = c(1.5,0.5),cex=2)
          text((0.75+0.35/2)*incubation.time,0.9,labels = paste0(round(min(variant.concentrations)/1e3)),adj = c(-1,0.5),cex=2)
        }
        if (scale.data==F){
          fields::colorbar.plot(x = 0.75*incubation.time,y = 0.9*diff(range(na.omit(c(data))))+min(na.omit(c(data))),strip = seq(0.3,0.9,0.6/11),col = grey.colors(12),strip.width = 0.05,strip.length = 0.35)
          text((0.75-0.35/2)*incubation.time,1*diff(range(na.omit(c(data))))+min(na.omit(c(data))),labels = paste0('[Decoy] (µM):'),adj = c(0.25,0.5),cex=2)
          text((0.75-0.35/2)*incubation.time,0.9*diff(range(na.omit(c(data))))+min(na.omit(c(data))),labels = paste0(round(max(variant.concentrations)/1e3)),adj = c(1.5,0.5),cex=2)
          text((0.75+0.35/2)*incubation.time,0.9*diff(range(na.omit(c(data))))+min(na.omit(c(data))),labels = paste0(round(min(variant.concentrations)/1e3)),adj = c(-1,0.5),cex=2)
        }
        title(paste0(enzyme,':  ',prey.molecule,' --> ',decoy.molecule),line = -5,outer = TRUE,cex.main = 4.5)

        # Plot decoy-dependence curves
        if (manual.fEP.adjust==T){
          cmp.models.coefficients.Mn.1=(cmp.models.coefficients.Mn.1*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),])))-FP.baseline)/(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-FP.baseline)
          cmp.models.coefficients.Mn.2=(cmp.models.coefficients.Mn.2*(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))+mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),])))-FP.baseline)/(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-FP.baseline)
        }
        cmp.models.coefficients.S=(1 - cmp.models.coefficients.Mn.1)/(2 - cmp.models.coefficients.Mn.2 - cmp.models.coefficients.Mn.1)
        cmp.models.coefficients.MnAVG=(cmp.models.coefficients.Mn.1+cmp.models.coefficients.Mn.2)/2
        S.true=rowMeans(cmp.models.coefficients.S,na.rm = TRUE)[1]
          cmp.models.coefficients.Mn.1x=(cmp.models.coefficients.Mn.1 + 2*S.true - 1)/(2*S.true)
          cmp.models.coefficients.Mn.2x=(-cmp.models.coefficients.Mn.2 + 2*S.true - 1)/(2*S.true - 2)
        k.off.1=cmp.models.coefficients.lambda.1*(1-cmp.models.coefficients.Mn.1x); if (coerce.offrates==T){k.off.1[k.off.1<=0]=0}
        k.off.2=cmp.models.coefficients.lambda.2*(1-cmp.models.coefficients.Mn.2x); if (coerce.offrates==T){k.off.2[k.off.2<=0]=0}
        kt.1=min(na.omit(k.off.1))
        kt.2=min(na.omit(k.off.2))
        par(fig=c(0.5,1,0.5,0.9),mar=c(4.5,5,3,1),oma = c(0,0,0,0),new=TRUE)
        plot(variant.concentrations/1e3,rowMeans(cmp.models.coefficients.S,na.rm = TRUE),main = 'Decoy vs Fast Events',xlab = paste0('[Decoy] (µM)'),ylab = expression('Proportion of k'[off]*''^'fast'), ylim = 0:1,type='p',cex.main=2,cex.lab=2,cex.axis=2)
        arrows(x0 = variant.concentrations/1e3, x1 = variant.concentrations/1e3, y0 = rowMeans(cmp.models.coefficients.S,na.rm = TRUE) - apply(cmp.models.coefficients.S,MARGIN = c(1),FUN = sd,na.rm = TRUE), y1 = rowMeans(cmp.models.coefficients.S,na.rm = TRUE) + apply(cmp.models.coefficients.S,MARGIN = c(1),FUN = sd,na.rm = TRUE),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        par(fig=c(0,0.35,0.25,0.5),mar=c(5,6,3,1),oma = c(0,0,0,0),new=TRUE)
        plot(variant.concentrations/1e3,rowMeans(cmp.models.coefficients.Mn.1x,na.rm = TRUE),main = expression('[EP]'['eq1']),xlab = paste0('[Decoy] (µM)'),ylab = expression('Fraction [EP]'[0]),type='p',cex.main=2,cex.lab=2,cex.axis=2)
        arrows(x0 = variant.concentrations/1e3, x1 = variant.concentrations/1e3, y0 = rowMeans(cmp.models.coefficients.Mn.1x,na.rm = TRUE) - apply(cmp.models.coefficients.Mn.1x,MARGIN = c(1),FUN = sd,na.rm = TRUE), y1 = rowMeans(cmp.models.coefficients.Mn.1x,na.rm = TRUE) + apply(cmp.models.coefficients.Mn.1x,MARGIN = c(1),FUN = sd,na.rm = TRUE),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        if (fit.DT.mdls==T){
          IC50a.data=list('y'=c(cmp.models.coefficients.Mn.1x),'D'=rep(variant.concentrations,times=ncol(cmp.models.coefficients.Mn.1x)))
          IC50a.mdl=IC50.fit(IC50a.data)
          try(lines(seq(0,max(variant.concentrations)*1e-3,1e-3),predict(IC50a.mdl,newdata=list('D' = seq(0,max(variant.concentrations),1))),lwd=3))
          if (show.constants==T & is.null(IC50a.mdl)==F){
            cc.temp = signif(coefficients(IC50a.mdl)[['IC50']],2)
            text(max(variant.concentrations)*1e-3*0.4,0.85,labels = substitute(paste('IC'[50],' = ',cc.temp,' nM',sep = "")),adj = c(0,0.5),cex=1)
            cc.temp = signif(coefficients(IC50a.mdl)[['h']],2)
            text(max(variant.concentrations)*1e-3*0.4,0.6,labels = substitute(paste('h = ',cc.temp,sep = "")),adj = c(0,0.5),cex=1)
          }
        }
        par(fig=c(0,0.35,0,0.25),mar=c(5,6,3,1),oma = c(0,0,0,0),new=TRUE)
        plot(variant.concentrations/1e3,rowMeans(cmp.models.coefficients.Mn.2x,na.rm = TRUE),main = expression('[EP]'['eq2']),xlab = paste0('[Decoy] (µM)'),ylab = expression('Fraction [EP]'[0]),type='p',cex.main=2,cex.lab=2,cex.axis=2)
        arrows(x0 = variant.concentrations/1e3, x1 = variant.concentrations/1e3, y0 = rowMeans(cmp.models.coefficients.Mn.2x,na.rm = TRUE) - apply(cmp.models.coefficients.Mn.2x,MARGIN = c(1),FUN = sd,na.rm = TRUE), y1 = rowMeans(cmp.models.coefficients.Mn.2x,na.rm = TRUE) + apply(cmp.models.coefficients.Mn.2x,MARGIN = c(1),FUN = sd,na.rm = TRUE),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        if (fit.DT.mdls==T){
          IC50b.data=list('y'=c(cmp.models.coefficients.Mn.2x),'D'=rep(variant.concentrations,times=ncol(cmp.models.coefficients.Mn.2x)))
          IC50b.mdl=IC50.fit(IC50b.data)
          try(lines(seq(0,max(variant.concentrations)*1e-3,1e-3),predict(IC50b.mdl,newdata=list('D' = seq(0,max(variant.concentrations),1))),lwd=3))
          if (show.constants==T & is.null(IC50b.mdl)==F){
            cc.temp = signif(coefficients(IC50b.mdl)[['IC50']],2)
            text(max(variant.concentrations)*1e-3*0.4,0.85,labels = substitute(paste('IC'[50],' = ',cc.temp,' nM',sep = "")),adj = c(0,0.5),cex=1)
            cc.temp = signif(coefficients(IC50b.mdl)[['h']],2)
            text(max(variant.concentrations)*1e-3*0.4,0.6,labels = substitute(paste('h = ',cc.temp,sep = "")),adj = c(0,0.5),cex=1)
          }
        }
        par(fig=c(0.35,1,0.25,0.5),mar=c(5,6,3,1),oma = c(0,0,0,0),new=TRUE)
        plot(variant.concentrations*1e-3,rowMeans(k.off.1,na.rm = TRUE)*1e3,main = expression('Decoy-Dependence of k'[off]*''^obs*''[fast]),xlab = paste0('[Decoy] (µM)'),ylab = expression('k'['off']*''^obs*''[fast]*' x10'^-3*' (s'^-1*')'),type='p',ylim=c(min(rowMeans(k.off.1,na.rm = TRUE)*1e3,na.rm = TRUE),1.5*as.numeric(rev(rowMeans(k.off.1,na.rm = TRUE)[order(rowMeans(k.off.1,na.rm = TRUE))])[2])*1e3),cex.main=3,cex.lab=2.5,cex.axis=2.5)
        arrows(x0 = variant.concentrations*1e-3, x1 = variant.concentrations*1e-3, y0 = rowMeans(k.off.1,na.rm = TRUE)*1e3 - apply(k.off.1,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3, y1 = rowMeans(k.off.1,na.rm = TRUE)*1e3 + apply(k.off.1,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3,code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        if (fit.DT.mdls==T){
          # Generate competition models
          DTmod1.data=list('x'=(rep(variant.concentrations,times=dim(data)[3]))[is.na(c(k.off.1-kt.1))==F],'y'=na.omit(c(k.off.1-kt.1)))
          DTmod1.both=DT.fit(DTmod1.data,match.optimizers)
          DTmod1.opt=DTmod1.both[['NULL']]
          DTmod1.opt2=DTmod1.both[['DT']]
          lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod1.opt,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt.1*1e3,col='green',lwd=3)
          lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod1.opt2,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt.1*1e3,col='purple',lwd=3)
          legend(legend.location,legend = c('Classic Competition Model','Direct Transfer Model'),col = c('green','purple'),fill = c('green','purple'),cex = 1)
        }
        par(fig=c(0.35,1,0,0.25),mar=c(5,6,3,1),oma = c(0,0,0,0),new=TRUE)
        plot(variant.concentrations*1e-3,rowMeans(k.off.2,na.rm = TRUE)*1e3,main = expression('Decoy-Dependence of k'[off]*''^obs*''[slow]),xlab = paste0('[Decoy] (µM)'),ylab = expression('k'['off']*''^obs*''[slow]*' x10'^-3*' (s'^-1*')'),type='p',ylim=c(min(rowMeans(k.off.2,na.rm = TRUE)*1e3,na.rm = TRUE),1.5*as.numeric(rev(rowMeans(k.off.2,na.rm = TRUE)[order(rowMeans(k.off.2,na.rm = TRUE))])[2])*1e3),cex.main=3,cex.lab=2.5,cex.axis=2.5)
        arrows(x0 = variant.concentrations*1e-3, x1 = variant.concentrations*1e-3, y0 = rowMeans(k.off.2,na.rm = TRUE)*1e3 - apply(k.off.2,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3, y1 = rowMeans(k.off.2,na.rm = TRUE)*1e3 + apply(k.off.2,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3,code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        if (fit.DT.mdls==T){
          # Generate competition models
          DTmod2.data=list('x'=(rep(variant.concentrations,times=dim(data)[3]))[is.na(c(k.off.2-kt.2))==F],'y'=na.omit(c(k.off.2-kt.2)))
          DTmod2.both=DT.fit(DTmod2.data,match.optimizers)
          DTmod2.opt=DTmod2.both[['NULL']]
          DTmod2.opt2=DTmod2.both[['DT']]
          lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod2.opt,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt.2*1e3,col='green',lwd=3)
          lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod2.opt2,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt.2*1e3,col='purple',lwd=3)
          legend(legend.location,legend = c('Classic Competition Model','Direct Transfer Model'),col = c('green','purple'),fill = c('green','purple'),cex = 1)
        }
        if (fit.DT.mdls==T){
          par(mfrow=c(2,1),mar=c(5,6,3,1),oma = c(0,0,0,0))
          plot(variant.concentrations*1e-3,rowMeans(k.off.1,na.rm = TRUE)*1e3,main = expression('Decoy-Dependence of k'[off]*''^obs*''[fast]),xlab = paste0('[Decoy] (µM)'),ylab = expression('k'['off']*''^obs*''[fast]*' x10'^-3*' (s'^-1*')'),type='p',ylim=c(min(rowMeans(k.off.1,na.rm = TRUE)*1e3,na.rm = TRUE),1.5*as.numeric(rev(rowMeans(k.off.1,na.rm = TRUE)[order(rowMeans(k.off.1,na.rm = TRUE))])[2])*1e3),cex.main=3,cex.lab=2.5,cex.axis=2.5)
          arrows(x0 = variant.concentrations*1e-3, x1 = variant.concentrations*1e-3, y0 = rowMeans(k.off.1,na.rm = TRUE)*1e3 - apply(k.off.1,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3, y1 = rowMeans(k.off.1,na.rm = TRUE)*1e3 + apply(k.off.1,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3,code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
          lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod1.opt,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt.1*1e3,col='green',lwd=3)
          lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod1.opt2,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt.1*1e3,col='purple',lwd=3)
          legend(legend.location,legend = c('Classic Competition Model','Direct Transfer Model'),col = c('green','purple'),fill = c('green','purple'),cex = 1.6)
          if (show.constants==T){
            cc.temp = signif((summary(DTmod1.opt)[['coefficients']][1,1])*1e3,2); text(0,(diff(range(rowMeans(k.off.1,na.rm = TRUE)*1e3))*0.5+min(rowMeans(k.off.1,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=1,col='green')
            cc.temp = signif((summary(DTmod1.opt2)[['coefficients']][1,1])*1e3,2); text(0,(diff(range(rowMeans(k.off.1,na.rm = TRUE)*1e3))*0.4+min(rowMeans(k.off.1,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=1,col='purple')
            cc.temp = signif((summary(DTmod1.opt2)[['coefficients']][4,1]),2); text(0,(diff(range(rowMeans(k.off.1,na.rm = TRUE)*1e3))*0.3+min(rowMeans(k.off.1,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[theta],' = ',cc.temp,' M'^-1,'s'^-1,sep = "")),adj = c(0,0.5),cex=1,col='purple')
          }
          plot(variant.concentrations*1e-3,rowMeans(k.off.2,na.rm = TRUE)*1e3,main = expression('Decoy-Dependence of k'[off]*''^obs*''[slow]),xlab = paste0('[Decoy] (µM)'),ylab = expression('k'['off']*''^obs*''[slow]*' x10'^-3*' (s'^-1*')'),type='p',ylim=c(min(rowMeans(k.off.2,na.rm = TRUE)*1e3,na.rm = TRUE),1.5*as.numeric(rev(rowMeans(k.off.2,na.rm = TRUE)[order(rowMeans(k.off.2,na.rm = TRUE))])[2])*1e3),cex.main=3,cex.lab=2.5,cex.axis=2.5)
          arrows(x0 = variant.concentrations*1e-3, x1 = variant.concentrations*1e-3, y0 = rowMeans(k.off.2,na.rm = TRUE)*1e3 - apply(k.off.2,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3, y1 = rowMeans(k.off.2,na.rm = TRUE)*1e3 + apply(k.off.2,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3,code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
          lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod2.opt,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt.2*1e3,col='green',lwd=3)
          lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(DTmod2.opt2,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt.2*1e3,col='purple',lwd=3)
          legend(legend.location,legend = c('Classic Competition Model','Direct Transfer Model'),col = c('green','purple'),fill = c('green','purple'),cex = 1.6)
          if (show.constants==T){
            cc.temp = signif((summary(DTmod2.opt)[['coefficients']][1,1])*1e3,2); text(0,(diff(range(rowMeans(k.off.2,na.rm = TRUE)*1e3))*0.5+min(rowMeans(k.off.2,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=1,col='green')
            cc.temp = signif((summary(DTmod2.opt2)[['coefficients']][1,1])*1e3,2); text(0,(diff(range(rowMeans(k.off.2,na.rm = TRUE)*1e3))*0.4+min(rowMeans(k.off.2,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=1,col='purple')
            cc.temp = signif((summary(DTmod2.opt2)[['coefficients']][4,1]),2); text(0,(diff(range(rowMeans(k.off.2,na.rm = TRUE)*1e3))*0.3+min(rowMeans(k.off.2,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[theta],' = ',cc.temp,' M'^-1,'s'^-1,sep = "")),adj = c(0,0.5),cex=1,col='purple')
          }
        }

        if (fit.DT.mdls==T){
          # Report and compare models
          if (is.null(IC50a.mdl)==F){
            print(paste('Equilibrium Competition 1 Summary:',sep = ''),quote = F); try(print(summary(IC50a.mdl),quote = F)); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          if (is.null(IC50b.mdl)==F){
            print(paste('Equilibrium Competition 2 Summary:',sep = ''),quote = F); try(print(summary(IC50b.mdl),quote = F)); print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          print(paste('(FAST) Classic Model Summary:',sep = ''),quote = F); print(summary(DTmod1.opt),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          print(paste('(FAST) Direct Transfer Model Summary:',sep = ''),quote = F); print(summary(DTmod1.opt2),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          fast.comparison=BIC(DTmod1.opt,DTmod1.opt2)
          print(paste("(FAST) Models' BIC Statistics:",sep = ''),quote = F); print(fast.comparison,quote = F)
          fast.BIC=BIC(DTmod1.opt2)-BIC(DTmod1.opt)
          print(paste('(FAST) ΔBIC = ',fast.BIC,sep = ''),quote = F)
          print(paste(rep('#',times=100),collapse = ''),quote = F)
          if (fast.BIC<=-10){
            print('These data favor the DIRECT TRANSFER (FAST) model!',quote = F)
            print(paste0('k.n1 (FAST) = ',signif(summary(DTmod1.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod1.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print(paste0('k.theta (FAST) = ',signif(summary(DTmod1.opt2)[['coefficients']][2,1],2),' ± ',signif(summary(DTmod1.opt2)[['coefficients']][2,2],2),' 1/M/s'),quote=F)
            print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          if (fast.BIC>=5){
            print('These data favor the CLASSIC (FAST) model!',quote = F)
            print(paste0('k.n1 (FAST) = ',signif(summary(DTmod1.opt)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod1.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          if (fast.BIC>-10 & fast.BIC<5){
            print('EITHER (FAST) model may have produced these data!',quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
            print(paste0('k.n1 (FAST) = ',signif(summary(DTmod1.opt)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod1.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print('OR',quote=F)
            print(paste0('k.n1 (FAST) = ',signif(summary(DTmod1.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod1.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print(paste0('k.theta (FAST) = ',signif(summary(DTmod1.opt2)[['coefficients']][2,1],2),' ± ',signif(summary(DTmod1.opt2)[['coefficients']][2,2],2),' 1/M/s'),quote=F)
            print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          #
          print(paste('(SLOW) Classic Model Summary:',sep = ''),quote = F); print(summary(DTmod2.opt),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          print(paste('(SLOW) Direct Transfer Model Summary:',sep = ''),quote = F); print(summary(DTmod2.opt2),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          slow.comparison=BIC(DTmod2.opt,DTmod2.opt2)
          print(paste("(SLOW) Models' BIC Statistics:",sep = ''),quote = F); print(slow.comparison,quote = F)
          slow.BIC=BIC(DTmod2.opt2)-BIC(DTmod2.opt)
          print(paste('(SLOW) ΔBIC = ',slow.BIC,sep = ''),quote = F)
          print(paste(rep('#',times=100),collapse = ''),quote = F)
          if (slow.BIC<=-10){
            print('These data favor the DIRECT TRANSFER (SLOW) model!',quote = F)
            print(paste0('k.n1 (SLOW) = ',signif(summary(DTmod2.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod2.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print(paste0('k.theta (SLOW) = ',signif(summary(DTmod2.opt2)[['coefficients']][2,1],2),' ± ',signif(summary(DTmod2.opt2)[['coefficients']][2,2],2),' 1/M/s'),quote=F)
            print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          if (slow.BIC>=5){
            print('These data favor the CLASSIC (SLOW) model!',quote = F)
            print(paste0('k.n1 (SLOW) = ',signif(summary(DTmod2.opt)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod2.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
          if (slow.BIC>-10 & slow.BIC<5){
            print('EITHER (SLOW) model may have produced these data!',quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
            print(paste0('k.n1 (SLOW) = ',signif(summary(DTmod2.opt)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod2.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print('OR',quote=F)
            print(paste0('k.n1 (SLOW) = ',signif(summary(DTmod2.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(DTmod2.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
            print(paste0('k.theta (SLOW) = ',signif(summary(DTmod2.opt2)[['coefficients']][2,1],2),' ± ',signif(summary(DTmod2.opt2)[['coefficients']][2,2],2),' 1/M/s'),quote=F)
            print(paste(rep('#',times=100),collapse = ''),quote = F)
          }
        }
      }
    }
  }
  if (experiment.type=='STOICH'){
    nlm.Kd=coefficients(BIND.models[['STOICH']])$Kd
    nlm.S=coefficients(BIND.models[['STOICH']])$S
    nlm.Mn=coefficients(BIND.models[['STOICH']])$Mn
    nlm.Mx=coefficients(BIND.models[['STOICH']])$Mx
    show(summary(BIND.models[['STOICH']]))

    scale.dat.shi=(data.scaled.eq-nlm.Mn)/(nlm.Mx-nlm.Mn)
    nlm.sims=data.frame('y'=predict(BIND.models[['STOICH']],newdata = list('E'=min(variant.concentrations):max(variant.concentrations))),'x'=min(variant.concentrations):max(variant.concentrations))
    nlm.sims.scaled=nlm.sims; nlm.sims.scaled$y=(nlm.sims$y-nlm.Mn)/(nlm.Mx-nlm.Mn)
    nlm.int=total.P/nlm.S

    plot(NULL,NULL,xlim=range(variant.concentrations),ylim=c(0,1),xlab=paste0('[',enzyme,'] (nM)'),ylab='Fraction Bound',main=paste(enzyme,' ',prey.molecule,'-Binding Stoichiometry'),cex.lab=1.5,cex.main=2,cex.axis=1.5)
    points(variant.concentrations,rowMeans(scale.dat.shi))
    lines(nlm.sims.scaled$x,nlm.sims.scaled$y)
    arrows(x0 = variant.concentrations,x1 = variant.concentrations,y0 = rowMeans(scale.dat.shi)-apply(scale.dat.shi,1,sd), y1 = rowMeans(scale.dat.shi)+apply(scale.dat.shi,1,sd),code = 3,lwd = 1,angle = 90,length = 0.5)

    abline(a=0,b=1/nlm.int,lty='solid',col='grey',lwd=3)
    abline(v=total.P,lty='dashed',col='blue',lwd=2)
    abline(v=nlm.int,lty='dashed',col='green',lwd=2)
    legend('bottomright',legend=c('[Ligand]','[Protein]'),col = c('blue','green'),fill = c('blue','green'),cex=2)
    text(x=0.5*max(variant.concentrations),y=0.4,labels=paste0('Ligand:Protein Stoichiometry = ',signif(nlm.S,2)),adj=c(0,1),col='red',cex=1.8)

  }

  if (plot.pdf==T){
    dev.off()
  }
  if (save.console==T){
    sink()
  }

  ### Save Results

  RESULTS=mget(ls())
  if (save.data==T){
    save(RESULTS,file = paste(path.to.file,save.name,'.RData',sep = ''))
  }
  return(RESULTS)

  ##### END SCRIPT

}

file.parse <- function(data.rows,background.rows=NULL,exp.type='full',par.file='par.csv',perp.file='perp.csv',path.to.file='./'){

  setwd(path.to.file)
  raw.par=read.csv(file = par.file,sep = ',',header = TRUE)
  raw.perp=read.csv(file = perp.file,sep = ',',header = TRUE)
  par.set=list(NULL)
  perp.set=list(NULL)
  if (exp.type=='single'){
    exp.n=length(data.rows)*2
    for (i in 1:exp.n){
      if (is.null(background.rows)==F){
        back.set=paste0(background.rows[i],1:24)
      }
      if (i%%2==1){
        data.set=paste0(data.rows[round(i/2+0.01)],1:12)
      }
      if (i%%2==0){
        data.set=paste0(data.rows[round(i/2+0.01)],13:24)
      }
      whole.set=c('Time..s.',data.set,if(is.null(background.rows)==F){back.set})
      par.set[[i]]=raw.par[,(colnames(raw.par)%in%whole.set)]
      perp.set[[i]]=raw.perp[,(colnames(raw.perp)%in%whole.set)]
      write.table(x = par.set[[i]],file = paste0('par',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
      write.table(x = perp.set[[i]],file = paste0('perp',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
    }
  }
  if (exp.type=='half'){
    exp.n=length(data.rows)
    for (i in 1:exp.n){
      if (is.null(background.rows)==F){
        back.set=paste0(background.rows[i],1:24)
      }
      data.set=paste0(data.rows[i],1:24)
      whole.set=c('Time..s.',data.set,if(is.null(background.rows)==F){back.set})
      par.set[[i]]=raw.par[,(colnames(raw.par)%in%whole.set)]
      perp.set[[i]]=raw.perp[,(colnames(raw.perp)%in%whole.set)]
      write.table(x = par.set[[i]],file = paste0('par',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
      write.table(x = perp.set[[i]],file = paste0('perp',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
    }
  }
  if (exp.type=='full'){
    exp.n=length(data.rows)/2
    for (i in 1:exp.n){
      if (is.null(background.rows)==F){
        back.set=paste0(background.rows[i],1:24)
      }
      data.set=c(paste0(data.rows[2*(i-1)+1],1:24),paste0(data.rows[2*(i-1)+2],1:24))
      whole.set=c('Time..s.',data.set,if(is.null(background.rows)==F){back.set})
      par.set[[i]]=raw.par[,(colnames(raw.par)%in%whole.set)]
      perp.set[[i]]=raw.perp[,(colnames(raw.perp)%in%whole.set)]
      write.table(x = par.set[[i]],file = paste0('par',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
      write.table(x = perp.set[[i]],file = paste0('perp',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
    }
  }

}

FPmultalyze <- function(experiment.type='Kd',path.to.file='./',file.range=NULL,file.name=NULL,enzyme=NULL,prey.molecule=NULL,decoy.molecule=NULL,parameters.file=NULL,save.data=T,save.name=NULL,plot.pdf=T,plot.name=NULL,save.console=T,save.id=NULL,...){

  if (is.null(file.range)==T){
    file.n=length(list.files(path = path.to.file,pattern = 'perp'))-1
    file.range=1:file.n
  }
  if (is.null(file.name)==T){
    file.name=cbind(paste0('par',file.range,'.txt'),paste0('perp',file.range,'.txt'))
  }
  if (is.null(enzyme)==T){
    enzyme=paste0('Predator',file.range)
  }
  if (is.null(prey.molecule)==T){
    prey.molecule=paste0('Prey',file.range)
  }
  if (is.null(decoy.molecule)==T){
    decoy.molecule=paste0('Decoy',file.range)
  }
  if (is.null(save.name)==T){
    save.name=paste0('Results',file.range)
  }
  if (is.null(plot.name)==T){
    plot.name=paste0('Graphs',file.range)
  }
  if (is.null(save.id)==T){
    save.id=paste0('Console',file.range)
  }
  for (i in 1:length(file.range)){
    FPalyze(experiment.type,path.to.file,file.name[i,],enzyme[i],prey.molecule[i],decoy.molecule[i],parameters.file,save.data,save.name[i],plot.pdf,plot.name[i],save.console,save.id[i],...)
  }

}

