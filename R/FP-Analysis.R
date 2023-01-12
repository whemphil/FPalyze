# FP_Analysis.R

# A script to analyze florescence polarization data for protein-nucleic acid binding and competition

###############

##### DEVELOPER NOTES

# None

###############

FPalyze <- function(experiment.type,path.to.file='./',file.name=c('par.txt','perp.txt'),enzyme='Predator',prey.molecule='Prey',decoy.molecule='Decoy',parameters.file=NULL,save.data=T,save.name='Results',plot.pdf=F,plot.name='Graphs',save.console=T,save.id='Console',
                    variant.concentrations=NULL,
                    data.size='full',
                    time.step=30,
                    t.zero=1.5,
                    incubation.time=NULL,
                    coerce.timepoints=F,
                    equilibrium.points=10,
                    outliers='none',
                    scale.data=F,
                    default.mP.values=F,
                    background.subtraction=F,
                    G.factor=1,
                    estimate.initials=T,
                    reverse.timepoints=F,
                    regression.approach='both',
                    estimated.kd=NULL,
                    regress.data=T,
                    coerce.offrates=T,
                    start.kn1=NULL,
                    start.ktheta=NULL,
                    start.Bhalf=NULL,
                    start.N=NULL,
                    match.optimizers=T,
                    manual.fEP.adjust=F,
                    FP.baseline=NULL,
                    legend.location='bottomright',
                    total.P=NULL,
                    estimated.S=NULL,
                    estimated.Mn=NULL,
                    estimated.Mx=NULL,
                    show.constants=T,
                    use.anisotropy=T){

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
    #
    if ((experiment.type=='Kd' | experiment.type=='STOICH') & is.null(incubation.time)==T){
      incubation.time=30 # window of data collection, min
    }
    if (experiment.type=='COMP' & is.null(incubation.time)==T){
      incubation.time=120 # window of data collection, min
    }
    #
    if (experiment.type=='Kd' & is.null(total.P)==T & regression.approach=='quad'){
      total.P=5 # total concentration of prey molecule, nM
    }
    if (experiment.type=='STOICH' & is.null(total.P)==T){
      total.P=250 # total concentration of prey molecule, nM
    }

  }

  ### Data Import and Analysis

  # Generate time vector
  t=seq(t.zero*60,t.zero*60+incubation.time*60,time.step)

  # Import and clean raw data
  if (default.mP.values==T){
    if (data.size=='full'){
      data=array(0,dim = c(length(t),12,4))
    }
    if (data.size=='half'){
      data=array(0,dim = c(length(t),12,2))
    }
    if (data.size=='single'){
      data=array(0,dim = c(length(t),12,1))
    }
    if (coerce.timepoints==F){
      raw=as.data.frame(read.csv(file = paste(path.to.file,file.name,sep = ''),header = TRUE,sep = '\t'))
    }
    if (coerce.timepoints==T){
      raw=as.data.frame(read.csv(file = paste(path.to.file,file.name,sep = ''),header = TRUE,sep = '\t'))[1:length(t),]
    }
    if (outliers[1]!='none'){
        raw[,colnames(raw) %in% outliers]=NA
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
  }
  if (default.mP.values==F){
    if (coerce.timepoints==F){
      raw.par=as.data.frame(read.csv(file = paste(path.to.file,file.name[1],sep = ''),header = TRUE,sep = '\t'))
      raw.perp=as.data.frame(read.csv(file = paste(path.to.file,file.name[2],sep = ''),header = TRUE,sep = '\t'))
    }
    if (coerce.timepoints==T){
      raw.par=as.data.frame(read.csv(file = paste(path.to.file,file.name[1],sep = ''),header = TRUE,sep = '\t'))[1:length(t),]
      raw.perp=as.data.frame(read.csv(file = paste(path.to.file,file.name[2],sep = ''),header = TRUE,sep = '\t'))[1:length(t),]
    }
    if (outliers[1]!='none'){
        raw.par[,colnames(raw.par) %in% outliers]=NA
        raw.perp[,colnames(raw.perp) %in% outliers]=NA
    }
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
    data.scaled=(data-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))/(mean(na.omit(c(data[1:3,which.min(variant.concentrations),])))-mean(na.omit(c(data[(length(t)-equilibrium.points+1):length(t),which.max(variant.concentrations),]))))
  }

  if (experiment.type=='Kd' | experiment.type=='STOICH'){
    # Calculate equilibrium values from scaled data
    if (reverse.timepoints==F & data.size!='single'){
      data.scaled.eq=apply(data.scaled[(length(t)-equilibrium.points+1):length(t),,],c(2,3),mean)
      if (regression.approach=='none'){
          data.eq=apply(data[(length(t)-equilibrium.points+1):length(t),,],c(2,3),mean)
      }
    }
    if (reverse.timepoints==T & data.size!='single'){
      data.scaled.eq=apply(data.scaled[1:equilibrium.points,,],c(2,3),mean)
      if (regression.approach=='none'){
          data.eq=apply(data.scaled[1:equilibrium.points,,],c(2,3),mean)
      }
    }
    if (reverse.timepoints==F & data.size=='single'){
      data.scaled.eq=apply(data.scaled[(length(t)-equilibrium.points+1):length(t),,],c(2),mean)
      if (regression.approach=='none'){
        data.eq=apply(data[(length(t)-equilibrium.points+1):length(t),,],c(2),mean)
      }
    }
    if (reverse.timepoints==T & data.size=='single'){
      data.scaled.eq=apply(data.scaled[1:equilibrium.points,,],c(2),mean)
      if (regression.approach=='none'){
        data.eq=apply(data.scaled[1:equilibrium.points,,],c(2),mean)
      }
    }

    # Fit equilibrium binding curve data
    if (experiment.type=='STOICH' | regression.approach=='quad'){
      if (data.size=='full'){model.data=list('E'=rep(variant.concentrations,times=4)[is.na(c(data.scaled.eq))==FALSE],'FP'=c(data.scaled.eq)[is.na(c(data.scaled.eq))==FALSE],'Pt'=total.P)}
      if (data.size=='half'){model.data=list('E'=rep(variant.concentrations,times=2)[is.na(c(data.scaled.eq))==FALSE],'FP'=c(data.scaled.eq)[is.na(c(data.scaled.eq))==FALSE],'Pt'=total.P)}
      if (data.size=='single'){model.data=list('E'=rep(variant.concentrations,times=1)[is.na(c(data.scaled.eq))==FALSE],'FP'=c(data.scaled.eq)[is.na(c(data.scaled.eq))==FALSE],'Pt'=total.P)}
    } else {
      if (data.size=='full'){model.data=list('E'=rep(variant.concentrations,times=4)[is.na(c(data.scaled.eq))==FALSE],'FP'=c(data.scaled.eq)[is.na(c(data.scaled.eq))==FALSE])}
      if (data.size=='half'){model.data=list('E'=rep(variant.concentrations,times=2)[is.na(c(data.scaled.eq))==FALSE],'FP'=c(data.scaled.eq)[is.na(c(data.scaled.eq))==FALSE])}
      if (data.size=='single'){model.data=list('E'=rep(variant.concentrations,times=1)[is.na(c(data.scaled.eq))==FALSE],'FP'=c(data.scaled.eq)[is.na(c(data.scaled.eq))==FALSE])}
    }
    if (estimate.initials==T & data.size!='single'){
      if (regression.approach!='quad' & experiment.type=='Kd'){
        estimated.kd=variant.concentrations[which.min(abs(rowMeans(data.scaled.eq)-(max(rowMeans(data.scaled.eq))-(max(rowMeans(data.scaled.eq))-min(rowMeans(data.scaled.eq)))/2)))]
      }
      if (regression.approach=='quad' | experiment.type=='STOICH'){
        estimated.kd=abs(variant.concentrations[which.min(abs(rowMeans(data.scaled.eq)-(max(rowMeans(data.scaled.eq))-(max(rowMeans(data.scaled.eq))-min(rowMeans(data.scaled.eq)))/2)))]-0.5*total.P)
        estimated.S=1
        estimated.Mn=min(model.data[['FP']])
        estimated.Mx=max(model.data[['FP']])
      }
    }
    if (estimate.initials==T & data.size=='single'){
      if (regression.approach!='quad' & experiment.type=='Kd'){
        estimated.kd=variant.concentrations[which.min(abs(data.scaled.eq-(max(data.scaled.eq)-(max(data.scaled.eq)-min(data.scaled.eq))/2)))]
      }
      if (regression.approach=='quad' | experiment.type=='STOICH'){
        estimated.kd=abs(variant.concentrations[which.min(abs(data.scaled.eq-(max(data.scaled.eq)-(max(data.scaled.eq)-min(data.scaled.eq))/2)))]-0.5*total.P)
        estimated.S=1
        estimated.Mn=min(model.data[['FP']])
        estimated.Mx=max(model.data[['FP']])
      }
    }
    if (regression.approach=='std'){
      model=nls(FP~(E/(E+Kd))*(Mx-Mn)+Mn,start = list(Kd=estimated.kd,Mx=max(model.data[['FP']]),Mn=min(model.data[['FP']])),data = model.data,control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE))
    }
    if (regression.approach=='hill'){
      model=nls(FP~(E^n/(E^n+Kd))*(Mx-Mn)+Mn,start = list(Kd=estimated.kd,Mx=max(model.data[['FP']]),Mn=min(model.data[['FP']]),n=1),data = model.data,control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE))
    }
    if (regression.approach=='both'){
      model.std=nls(FP~(E/(E+Kd))*(Mx-Mn)+Mn,start = list(Kd=estimated.kd,Mx=max(model.data[['FP']]),Mn=min(model.data[['FP']])),data = model.data,control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE))
      model.hill=nls(FP~(E^n/(E^n+Kd))*(Mx-Mn)+Mn,start = list(Kd=estimated.kd,Mx=max(model.data[['FP']]),Mn=min(model.data[['FP']]),n=1),data = model.data,control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE))
    }
    if (regression.approach=='quad'){
      model.quad=nls(FP~(Mx-Mn)*(((Pt+(E+Kd))-((Pt+(E+Kd))^2-4*Pt*E)^0.5)/(2*Pt))+Mn,model.data,list('Kd'=estimated.kd,'Mn'=min(model.data[['FP']]),'Mx'=max(model.data[['FP']])))
    }
    if (experiment.type=='STOICH'){
      model.quadS=nls(FP~(Mx-Mn)*(((Pt+S*(E+Kd))-((Pt+S*(E+Kd))^2-4*Pt*S*E)^0.5)/(2*Pt))+Mn,model.data,list('Kd'=estimated.kd,'S'=estimated.S,'Mn'=estimated.Mn,'Mx'=estimated.Mx))
    }
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

  if (experiment.type=='Kd' & data.size=='full'){
    # Plot association curves
    par(mfrow=c(4,3))
    for (i in 1:12){
      if (scale.data==T){
        plot(rep(t/60,times=4),c(data.scaled[,i,]),type = 'p',main = paste(enzyme,'-',prey.molecule,' Association Curve:  ',enzyme,' at ',variant.concentrations[i],' nM',sep = ''),xlab = 'Time (minutes)',ylab = 'Relative Polarization',ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex.main=0.8,cex=1,col='black')
      }
      if (scale.data==F){
        plot(rep(t/60,times=4),c(data.scaled[,i,]),type = 'p',main = paste(enzyme,'-',prey.molecule,' Association Curve:  ',enzyme,' at ',variant.concentrations[i],' nM',sep = ''),xlab = 'Time (minutes)',ylab = 'Polarization (mP)',ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex.main=0.8,cex=1,col='black')
      }
      lines(t/60,rowMeans(data.scaled[,i,],na.rm = TRUE),col='red',lwd=4)
    }

    if (regression.approach=='none'){
      # Plot equilibrium binding curve
      par(mfrow=c(1,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=4),c(data.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Unmodeled Binding Data',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization')
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=4),c(data.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Unmodeled Binding Data',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)')
      }
      lines(variant.concentrations,rowMeans(data.eq,na.rm = TRUE),col='blue',lty='dashed')
    }
    if (regression.approach=='std'){
      # Plot equilibrium binding curve
      par(mfrow=c(1,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=4),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Standard Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=4),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Standard Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      lines(min(variant.concentrations):max(variant.concentrations),predict(model,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))))

      # Summarize binding curve regression
      print(summary(model)); print(paste(rep('#',times=100),collapse = ''),quote = F)
    }
    if (regression.approach=='hill'){
      # Plot equilibrium binding curve
      par(mfrow=c(1,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=4),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Hill Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=4),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Hill Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      lines(min(variant.concentrations):max(variant.concentrations),predict(model,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))))

      # Summarize binding curve regression
      print(summary(model)); print(paste(rep('#',times=100),collapse = ''),quote = F)
    }
    if (regression.approach=='both'){
      # Plot equilibrium binding curve
      par(mfrow=c(2,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=4),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.std,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.hill,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.std,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.hill,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=4),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.std,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.hill,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.std,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.hill,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))
      }

      # Summarize binding curve regression
      print(summary(model.std)); print(paste(rep('#',times=100),collapse = ''),quote = F)
      print(summary(model.hill)); print(paste(rep('#',times=100),collapse = ''),quote = F)

      # Model comparison statistics
      BIC=BIC(model.std,model.hill)
      print(BIC,quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      delta.BIC=BIC(model.hill)-BIC(model.std)
      print(paste('ΔBIC = ',delta.BIC,sep = ''),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      if (delta.BIC<=-10){
        print('These data favor the HILL model!',quote = F)
        print(paste0('Kd = ',signif(summary(model.hill)[['coefficients']][1,1],2),' ± ',signif(summary(model.hill)[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste0('n = ',signif(summary(model.hill)[['coefficients']][4,1],2),' ± ',signif(summary(model.hill)[['coefficients']][4,2],2)),quote=F)
        print(paste0('B_0.5 ≈ ',signif(summary(model.hill)[['coefficients']][1,1]^(1/summary(model.hill)[['coefficients']][4,1]),2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (delta.BIC>=5){
        print('These data favor the STANDARD model!',quote = F)
        print(paste0('Kd = ',signif(summary(model.std)[['coefficients']][1,1],2),' ± ',signif(summary(model.std)[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (delta.BIC>-10 & delta.BIC<5){
        print('EITHER model may have produced these data!',quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        print(paste0('Kd = ',signif(summary(model.std)[['coefficients']][1,1],2),' ± ',signif(summary(model.std)[['coefficients']][1,2],2),' nM'),quote=F)
        print('OR',quote=F)
        print(paste0('Kd = ',signif(summary(model.hill)[['coefficients']][1,1],2),' ± ',signif(summary(model.hill)[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste0('n = ',signif(summary(model.hill)[['coefficients']][4,1],2),' ± ',signif(summary(model.hill)[['coefficients']][4,2],2)),quote=F)
        print(paste0('B_0.5 ≈ ',signif(summary(model.hill)[['coefficients']][1,1]^(1/summary(model.hill)[['coefficients']][4,1]),2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }

      # Outlier Identification
      if (outliers[1]=='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]=='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[-((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
    }
    if (regression.approach=='quad'){
      # Plot equilibrium binding curve
      par(mfrow=c(2,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=4),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.quad,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.quad,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=4),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.quad,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.quad,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
      }

      # Summarize binding curve regression
      print(summary(model.quad)); print(paste(rep('#',times=100),collapse = ''),quote = F)

      # Outlier Identification
      if (outliers[1]=='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]=='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[-((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
    }

  }
  if (experiment.type=='Kd' & data.size=='half'){
    # Plot association curves
    par(mfrow=c(4,3))
    for (i in 1:12){
      if (scale.data==T){
        plot(rep(t/60,times=2),c(data.scaled[,i,]),type = 'p',main = paste(enzyme,'-',prey.molecule,' Association Curve:  ',enzyme,' at ',variant.concentrations[i],' nM',sep = ''),xlab = 'Time (minutes)',ylab = 'Relative Polarization',ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex.main=0.8,cex=1,col='black')
      }
      if (scale.data==F){
        plot(rep(t/60,times=2),c(data.scaled[,i,]),type = 'p',main = paste(enzyme,'-',prey.molecule,' Association Curve:  ',enzyme,' at ',variant.concentrations[i],' nM',sep = ''),xlab = 'Time (minutes)',ylab = 'Polarization (mP)',ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex.main=0.8,cex=1,col='black')
      }
      lines(t/60,rowMeans(data.scaled[,i,],na.rm = TRUE),col='red',lwd=4)
    }

    if (regression.approach=='none'){
      # Plot equilibrium binding curve
      par(mfrow=c(1,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=2),c(data.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Unmodeled Binding Data',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization')
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=2),c(data.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Unmodeled Binding Data',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)')
      }
      lines(variant.concentrations,rowMeans(data.eq,na.rm = TRUE),col='blue',lty='dashed')
    }
    if (regression.approach=='std'){
      # Plot equilibrium binding curve
      par(mfrow=c(1,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Standard Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Standard Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      lines(min(variant.concentrations):max(variant.concentrations),predict(model,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))))

      # Summarize binding curve regression
      print(summary(model)); print(paste(rep('#',times=100),collapse = ''),quote = F)
    }
    if (regression.approach=='hill'){
      # Plot equilibrium binding curve
      par(mfrow=c(1,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Hill Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Hill Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      lines(min(variant.concentrations):max(variant.concentrations),predict(model,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))))

      # Summarize binding curve regression
      print(summary(model)); print(paste(rep('#',times=100),collapse = ''),quote = F)
    }
    if (regression.approach=='both'){
      # Plot equilibrium binding curve
      par(mfrow=c(2,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.std,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.hill,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.std,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.hill,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.std,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.hill,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.std,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.hill,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))
      }

      # Summarize binding curve regression
      print(summary(model.std)); print(paste(rep('#',times=100),collapse = ''),quote = F)
      print(summary(model.hill)); print(paste(rep('#',times=100),collapse = ''),quote = F)

      # Model comparison statistics
      BIC=BIC(model.std,model.hill)
      print(BIC,quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      delta.BIC=BIC(model.hill)-BIC(model.std)
      print(paste('ΔBIC = ',delta.BIC,sep = ''),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      if (delta.BIC<=-10){
        print('These data favor the HILL model!',quote = F)
        print(paste0('Kd = ',signif(summary(model.hill)[['coefficients']][1,1],2),' ± ',signif(summary(model.hill)[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste0('n = ',signif(summary(model.hill)[['coefficients']][4,1],2),' ± ',signif(summary(model.hill)[['coefficients']][4,2],2)),quote=F)
        print(paste0('B_0.5 ≈ ',signif(summary(model.hill)[['coefficients']][1,1]^(1/summary(model.hill)[['coefficients']][4,1]),2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (delta.BIC>=5){
        print('These data favor the STANDARD model!',quote = F)
        print(paste0('Kd = ',signif(summary(model.std)[['coefficients']][1,1],2),' ± ',signif(summary(model.std)[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (delta.BIC>-10 & delta.BIC<5){
        print('EITHER model may have produced these data!',quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        print(paste0('Kd = ',signif(summary(model.std)[['coefficients']][1,1],2),' ± ',signif(summary(model.std)[['coefficients']][1,2],2),' nM'),quote=F)
        print('OR',quote=F)
        print(paste0('Kd = ',signif(summary(model.hill)[['coefficients']][1,1],2),' ± ',signif(summary(model.hill)[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste0('n = ',signif(summary(model.hill)[['coefficients']][4,1],2),' ± ',signif(summary(model.hill)[['coefficients']][4,2],2)),quote=F)
        print(paste0('B_0.5 ≈ ',signif(summary(model.hill)[['coefficients']][1,1]^(1/summary(model.hill)[['coefficients']][4,1]),2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }

      # Outlier Identification
      if (outliers[1]=='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]=='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[-((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
    }
    if (regression.approach=='quad'){
      # Plot equilibrium binding curve
      par(mfrow=c(2,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.quad,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.quad,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.quad,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.quad,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
      }

      # Summarize binding curve regression
      print(summary(model.quad)); print(paste(rep('#',times=100),collapse = ''),quote = F)

      # Outlier Identification
      if (outliers[1]=='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]=='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[-((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
    }

  }

  if (experiment.type=='Kd' & data.size=='single'){
    # Plot association curves
    par(mfrow=c(4,3))
    for (i in 1:12){
      if (scale.data==T){
        plot(rep(t/60,times=1),c(data.scaled[,i,]),type = 'p',main = paste(enzyme,'-',prey.molecule,' Association Curve:  ',enzyme,' at ',variant.concentrations[i],' nM',sep = ''),xlab = 'Time (minutes)',ylab = 'Relative Polarization',ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex.main=0.8,cex=1,col='black')
      }
      if (scale.data==F){
        plot(rep(t/60,times=1),c(data.scaled[,i,]),type = 'p',main = paste(enzyme,'-',prey.molecule,' Association Curve:  ',enzyme,' at ',variant.concentrations[i],' nM',sep = ''),xlab = 'Time (minutes)',ylab = 'Polarization (mP)',ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex.main=0.8,cex=1,col='black')
      }
      lines(t/60,data.scaled[,i,],col='red',lwd=4)
    }

    if (regression.approach=='none'){
      # Plot equilibrium binding curve
      par(mfrow=c(1,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=1),c(data.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Unmodeled Binding Data',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization')
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=1),c(data.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Unmodeled Binding Data',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)')
      }
      lines(variant.concentrations,data.eq,col='blue',lty='dashed')
    }
    if (regression.approach=='std'){
      # Plot equilibrium binding curve
      par(mfrow=c(1,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=1),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Standard Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=1),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Standard Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      lines(min(variant.concentrations):max(variant.concentrations),predict(model,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))))

      # Summarize binding curve regression
      print(summary(model)); print(paste(rep('#',times=100),collapse = ''),quote = F)
    }
    if (regression.approach=='hill'){
      # Plot equilibrium binding curve
      par(mfrow=c(1,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=1),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Hill Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=1),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Hill Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(0,max(na.omit(c(data.scaled.eq)))))
      }
      lines(min(variant.concentrations):max(variant.concentrations),predict(model,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))))

      # Summarize binding curve regression
      print(summary(model)); print(paste(rep('#',times=100),collapse = ''),quote = F)
    }
    if (regression.approach=='both'){
      # Plot equilibrium binding curve
      par(mfrow=c(2,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=1),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.std,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.hill,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))

        plot(log10(variant.concentrations),data.scaled.eq,type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.std,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.hill,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=1),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.std,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.hill,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))

        plot(log10(variant.concentrations),data.scaled.eq,type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.std,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.hill,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='red',lty='dashed')
        legend('bottomright',legend = c('Standard Model','Hill Model'),fill = c('blue','red'),col = c('blue','red'))
      }

      # Summarize binding curve regression
      print(summary(model.std)); print(paste(rep('#',times=100),collapse = ''),quote = F)
      print(summary(model.hill)); print(paste(rep('#',times=100),collapse = ''),quote = F)

      # Model comparison statistics
      BIC=BIC(model.std,model.hill)
      print(BIC,quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      delta.BIC=BIC(model.hill)-BIC(model.std)
      print(paste('ΔBIC = ',delta.BIC,sep = ''),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      if (delta.BIC<=-10){
        print('These data favor the HILL model!',quote = F)
        print(paste0('Kd = ',signif(summary(model.hill)[['coefficients']][1,1],2),' ± ',signif(summary(model.hill)[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste0('n = ',signif(summary(model.hill)[['coefficients']][4,1],2),' ± ',signif(summary(model.hill)[['coefficients']][4,2],2)),quote=F)
        print(paste0('B_0.5 ≈ ',signif(summary(model.hill)[['coefficients']][1,1]^(1/summary(model.hill)[['coefficients']][4,1]),2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (delta.BIC>=5){
        print('These data favor the STANDARD model!',quote = F)
        print(paste0('Kd = ',signif(summary(model.std)[['coefficients']][1,1],2),' ± ',signif(summary(model.std)[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (delta.BIC>-10 & delta.BIC<5){
        print('EITHER model may have produced these data!',quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        print(paste0('Kd = ',signif(summary(model.std)[['coefficients']][1,1],2),' ± ',signif(summary(model.std)[['coefficients']][1,2],2),' nM'),quote=F)
        print('OR',quote=F)
        print(paste0('Kd = ',signif(summary(model.hill)[['coefficients']][1,1],2),' ± ',signif(summary(model.hill)[['coefficients']][1,2],2),' nM'),quote=F)
        print(paste0('n = ',signif(summary(model.hill)[['coefficients']][4,1],2),' ± ',signif(summary(model.hill)[['coefficients']][4,2],2)),quote=F)
        print(paste0('B_0.5 ≈ ',signif(summary(model.hill)[['coefficients']][1,1]^(1/summary(model.hill)[['coefficients']][4,1]),2),' nM'),quote=F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
      }

      # Outlier Identification
      if (outliers[1]=='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]=='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.hill[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[-((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
    }
    if (regression.approach=='quad'){
      # Plot equilibrium binding curve
      par(mfrow=c(2,1))
      if (scale.data==T){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.quad,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'Relative EQ-Polarization',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.quad,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
      }
      if (scale.data==F){
        plot(rep(variant.concentrations,times=2),c(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        lines(min(variant.concentrations):max(variant.concentrations),predict(model.quad,newdata=list('E'=min(variant.concentrations):max(variant.concentrations))),col='blue',lty='dashed')

        plot(log10(variant.concentrations),rowMeans(data.scaled.eq),type = 'p',main = paste(enzyme,'-',prey.molecule,' Binding Curve',sep = ''),xlab = paste('log10[',enzyme,'] (nM)',sep = ''),ylab = 'EQ-Polarization (mP)',ylim = c(min(na.omit(c(data.scaled.eq))),max(na.omit(c(data.scaled.eq)))))
        arrows(x0 = log10(variant.concentrations), x1 = log10(variant.concentrations), y0 = rowMeans(data.scaled.eq) - apply(data.scaled.eq,MARGIN = c(1),FUN = sd), y1 = rowMeans(data.scaled.eq) + apply(data.scaled.eq,MARGIN = c(1),FUN = sd),code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)
        lines(log10(seq(min(variant.concentrations),max(variant.concentrations),0.01)),predict(model.quad,newdata=list('E'=seq(min(variant.concentrations),max(variant.concentrations),0.01))),col='blue',lty='dashed')
      }

      # Summarize binding curve regression
      print(summary(model.quad)); print(paste(rep('#',times=100),collapse = ''),quote = F)

      # Outlier Identification
      if (outliers[1]=='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]=='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==F){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==T){
        temp.1=data.scaled.eq; temp.1[is.na(temp.1)==FALSE]=(environment(model.quad[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[-((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
    }

  }

  if (experiment.type=='COMP' & data.size=='full'){
    # Plot dissociation curves
    par(mfrow=c(4,3),mar=c(4.5,5,3,1),oma = c(0,0,5,0))
    cmp.models=list(NULL)
    cmp.models.coefficients.lambda=matrix(NA,ncol=4,nrow=length(variant.concentrations)); cmp.models.coefficients.Mn=matrix(NA,ncol=4,nrow=length(variant.concentrations))
    for (i in 1:12){
      if (scale.data==T){
        plot(rep(t/60,times=4),c(data.scaled[,i,]),type = 'p',main = paste0('[Decoy] = ',signif(variant.concentrations[i]/1e3,2),' µM'),xlab = 'Time (min)',ylab = 'Relative Polarization',ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex=1,cex.main=2.5,cex.axis=2,cex.lab=2)
      }
      if (scale.data==F){
        plot(rep(t/60,times=4),c(data[,i,]),type = 'p',main = paste0('[Decoy] = ',signif(variant.concentrations[i]/1e3,2),' µM'),xlab = 'Time (min)',ylab = 'Polarization (mP)',ylim = c(min(na.omit(c(data))),max(na.omit(c(data)))),cex=1,cex.main=2.5,cex.axis=2,cex.lab=2)
      }
      if (i==1){
        title(paste0(enzyme,':  ',prey.molecule,' --> ',decoy.molecule),line = 1.5,outer = TRUE,cex.main = 4)
      }
      if (regress.data==T){
        for (q in 1:4){
          cmp.models.data=list('FP'=c(data.scaled[,i,q]),'tt'=t); try(cmp.models[[paste(i,q,sep = '_')]]<-nls(FP~(1-Mn)*exp(-kn1*tt)+Mn,data=cmp.models.data,start=list('kn1'=1e-3,'Mn'=0),control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE)))
          if (is.null(cmp.models[[paste(i,q,sep = '_')]])==F){
            cmp.models.coefficients.lambda[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[1]; cmp.models.coefficients.Mn[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[2]
          }
          if (is.null(cmp.models[[paste(i,q,sep = '_')]])==T){
            cmp.models.coefficients.lambda[i,q]=0; cmp.models.coefficients.Mn[i,q]=1
            warning(paste0('Reaction ',i,'-',q,' does not fit an exponential decay -- it was coerced to a horizontal line at binding saturation.'))
          }
        }
        cmp.models.data=list('FP'=c(data.scaled[,i,]),'tt'=rep(t,times=4)); try(cmp.models[[paste(i,'all',sep = '_')]]<-nls(FP~(1-Mn)*exp(-kn1*tt)+Mn,data=cmp.models.data,start=list('kn1'=1e-3,'Mn'=0),control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE)))
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
        plot(NULL,NULL,main = 'Comparison of Competition Reactions',xlab = 'Time (min)',ylab = 'Relative Polarization',ylim = c(0,1),xlim = c(min(t)/60,max(t)/60),cex.main=2,cex.axis=2,cex.lab=2)
      }
      if (scale.data==F){
        plot(NULL,NULL,main = 'Comparison of Competition Reactions',xlab = 'Time (min)',ylab = 'Polarization (mP)',ylim = c(min(na.omit(c(data))),max(na.omit(c(data)))),xlim = c(min(t)/60,max(t)/60),cex.main=2,cex.axis=2,cex.lab=2)
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
      par(fig=c(0,1,0,0.5),mar=c(5,6,3,1),oma = c(0,0,0,0),new=TRUE)
      plot(variant.concentrations*1e-3,rowMeans(k.off,na.rm = TRUE)*1e3,main = expression('Decoy-Dependence of k'[off]*''^obs),xlab = paste0('[Decoy] (µM)'),ylab = expression('k'['off']*''^obs*' x10'^-3*' (s'^-1*')'),type='p',cex.main=3,cex.lab=2.5,cex.axis=2.5)
      arrows(x0 = variant.concentrations*1e-3, x1 = variant.concentrations*1e-3, y0 = rowMeans(k.off,na.rm = TRUE)*1e3 - apply(k.off,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3, y1 = rowMeans(k.off,na.rm = TRUE)*1e3 + apply(k.off,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3,code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)

      # Generate competition models
      fun.data=list('x'=(rep(variant.concentrations,times=4))[is.na(c(k.off-kt))==F],'y'=na.omit(c(k.off-kt)))
      if (estimate.initials==T){
        start.N=1
        start.ktheta=(rowMeans(k.off-kt,na.rm = TRUE)[1]-rowMeans(k.off-kt,na.rm = TRUE)[2])/(variant.concentrations[1]*1e-9-variant.concentrations[2]*1e-9)
        start.kn1=rowMeans(k.off-kt,na.rm = TRUE)[1]-(variant.concentrations[1]*1e-9)*start.ktheta
        start.Bhalf=variant.concentrations[which.min(abs(start.kn1/2-rowMeans(k.off-kt,na.rm = TRUE)))]
      }
      fun.model.opt2=nls(y~k.n1*(x*1e-9)^N/((x*1e-9)^N+(Bhalf*1e-9)^N)+k.theta*(x*1e-9),data = fun.data,start = list('k.n1'=start.kn1,'Bhalf'=start.Bhalf,'N'=start.N,'k.theta'=start.ktheta),control=list(minFactor=1e-10,maxiter=1e2,warnOnly=TRUE))
      if (match.optimizers==T & exists('fun.model.opt2')==T){
        fun.data<-list('x'=(rep(variant.concentrations,times=4))[is.na(c(k.off-kt))==F],'y'=na.omit(c(k.off-kt)),'Bhalf'=coefficients(fun.model.opt2)[['Bhalf']],'N'=coefficients(fun.model.opt2)[['N']])
        try(fun.model.opt<-nls(y~k.n1*(x*1e-9)^N/((x*1e-9)^N+(Bhalf*1e-9)^N),data = fun.data,start = list('k.n1'=coefficients(fun.model.opt2)[['k.n1']]),control=list(minFactor=1e-10,maxiter=1e2,warnOnly=TRUE)))
      }
      if (match.optimizers==F){
        fun.data=list('x'=(rep(variant.concentrations,times=4))[is.na(c(k.off-kt))==F],'y'=na.omit(c(k.off-kt)))
        try(fun.model.opt<-nls(y~k.n1*(x*1e-9)^N/((x*1e-9)^N+(Bhalf*1e-9)^N),data = fun.data,start = list('k.n1'=coefficients(fun.model.opt2)[['k.n1']],'Bhalf'=coefficients(fun.model.opt2)[['Bhalf']],'N'=coefficients(fun.model.opt2)[['N']]),control=list(minFactor=1e-10,maxiter=1e2,warnOnly=TRUE)))
      }
      if (exists('fun.model.opt')==T){
        lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(fun.model.opt,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt*1e3,col='green',lwd=3,lty='dashed')
        if (show.constants==T){
          cc.temp = signif((summary(fun.model.opt)[['coefficients']][1,1])*1e3,2); text(max(variant.concentrations)*1e-3*0.7,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.5+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=2.4,col='green')
        }
      }
      if (exists('fun.model.opt2')==T){
        lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(fun.model.opt2,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt*1e3,col='purple',lwd=3,lty='dashed')
        #legend(legend.location,legend = c('Classic Competition Model','Displacement-Transfer Model'),col = c('green','purple'),fill = c('green','purple'),cex = 3)
        legend(legend.location,legend = c('Classic Competition Model','Direct Transfer Model'),col = c('green','purple'),fill = c('green','purple'),cex = 3)
        if (show.constants==T){
          cc.temp = signif((summary(fun.model.opt2)[['coefficients']][1,1])*1e3,2); text(max(variant.concentrations)*1e-3*0.7,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.4+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=2.4,col='purple')
          cc.temp = signif((summary(fun.model.opt2)[['coefficients']][4,1]),2); text(max(variant.concentrations)*1e-3*0.7,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.3+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[theta],' = ',cc.temp,' M'^-1,'s'^-1,sep = "")),adj = c(0,0.5),cex=2.4,col='purple')
        }
      }

      # Report and compare models
      if (exists('fun.model.opt')==T){
        print(paste('Classic Model Summary:',sep = ''),quote = F); print(summary(fun.model.opt),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (exists('fun.model.opt2')==T){
        print(paste('Direct Transfer Model Summary:',sep = ''),quote = F); print(summary(fun.model.opt2),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (exists('fun.model.opt2')==T & exists('fun.model.opt')==T){
        model.comparison=BIC(fun.model.opt,fun.model.opt2)
        print(paste("Models' BIC Statistics:",sep = ''),quote = F); print(model.comparison,quote = F)
        delta.BIC=BIC(fun.model.opt2)-BIC(fun.model.opt)
        print(paste('ΔBIC = ',delta.BIC,sep = ''),quote = F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
        if (delta.BIC<=-10){
          print('These data favor the DIRECT TRANSFER model!',quote = F)
          print(paste0('k.n1 = ',signif(summary(fun.model.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(fun.model.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
          print(paste0('k.theta = ',signif(summary(fun.model.opt2)[['coefficients']][4,1],2),' ± ',signif(summary(fun.model.opt2)[['coefficients']][4,2],2),' 1/M/s'),quote=F)
          print(paste(rep('#',times=100),collapse = ''),quote = F)
          if (show.constants==T){
            text(max(variant.concentrations)*1e-3*0.65,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.05+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels='BIC Favors >',cex = 5,col = 'purple',adj=c(1.5,0.5))
          }
        }
        if (delta.BIC>=5){
          print('These data favor the CLASSIC model!',quote = F)
          print(paste0('k.n1 = ',signif(summary(fun.model.opt)[['coefficients']][1,1],2),' ± ',signif(summary(fun.model.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
          print(paste(rep('#',times=100),collapse = ''),quote = F)
          if (show.constants==T){
            text(max(variant.concentrations)*1e-3*0.65,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.15+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels='BIC Favors >',cex = 5,col = 'green',adj=c(1.5,0.5))
          }
        }
        if (delta.BIC>-10 & delta.BIC<5){
          print('EITHER model may have produced these data!',quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          print(paste0('k.n1 = ',signif(summary(fun.model.opt)[['coefficients']][1,1],2),' ± ',signif(summary(fun.model.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
          print('OR',quote=F)
          print(paste0('k.n1 = ',signif(summary(fun.model.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(fun.model.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
          print(paste0('k.theta = ',signif(summary(fun.model.opt2)[['coefficients']][4,1],2),' ± ',signif(summary(fun.model.opt2)[['coefficients']][4,2],2),' 1/M/s'),quote=F)
          print(paste(rep('#',times=100),collapse = ''),quote = F)
          if (show.constants==T){
            text(max(variant.concentrations)*1e-3*0.5,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.1+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels='BIC',cex = 5,col = 'grey',adj=c(1.5,0.5))
          }
        }
      }

      # Outlier Identification
      if (outliers[1]=='none' & default.mP.values==F){
        temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(fun.model.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]=='none' & default.mP.values==T){
        temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(fun.model.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==F){
        temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(fun.model.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[((colnames(raw.par)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==T){
        temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(fun.model.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))])[-((colnames(raw)[2:49])[c(seq(1,24,2),seq(2,24,2),seq(25,48,2),seq(26,48,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
    }
  }
  if (experiment.type=='COMP' & data.size=='half'){
    # Plot dissociation curves
    par(mfrow=c(4,3),mar=c(4.5,5,3,1),oma = c(0,0,5,0))
    cmp.models=list(NULL)
    cmp.models.coefficients.lambda=matrix(NA,ncol=2,nrow=length(variant.concentrations)); cmp.models.coefficients.Mn=matrix(NA,ncol=2,nrow=length(variant.concentrations))
    for (i in 1:12){
      if (scale.data==T){
        plot(rep(t/60,times=2),c(data.scaled[,i,]),type = 'p',main = paste0(decoy.molecule,' = ',signif(variant.concentrations[i]/1e3,2),' µM'),xlab = 'Time (min)',ylab = 'Relative Polarization',ylim = c(min(na.omit(c(data.scaled))),max(na.omit(c(data.scaled)))),cex=1,cex.main=2.5,cex.axis=2,cex.lab=2)
      }
      if (scale.data==F){
        plot(rep(t/60,times=2),c(data[,i,]),type = 'p',main = paste0(decoy.molecule,' = ',signif(variant.concentrations[i]/1e3,2),' µM'),xlab = 'Time (min)',ylab = 'Polarization (mP)',ylim = c(min(na.omit(c(data))),max(na.omit(c(data)))),cex=1,cex.main=2.5,cex.axis=2,cex.lab=2)
      }
      if (i==1){
        title(paste0(enzyme,':  ',prey.molecule,' --> ',decoy.molecule),line = 1.5,outer = TRUE,cex.main = 4)
      }
      if (regress.data==T){
        for (q in 1:2){
          cmp.models.data=list('FP'=c(data.scaled[,i,q]),'tt'=t); try(cmp.models[[paste(i,q,sep = '_')]]<-nls(FP~(1-Mn)*exp(-kn1*tt)+Mn,data=cmp.models.data,start=list('kn1'=1e-3,'Mn'=0),control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE)))
          if (is.null(cmp.models[[paste(i,q,sep = '_')]])==F){
            cmp.models.coefficients.lambda[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[1]; cmp.models.coefficients.Mn[i,q]=coefficients(cmp.models[[paste(i,q,sep = '_')]])[2]
          }
          if (is.null(cmp.models[[paste(i,q,sep = '_')]])==T){
            cmp.models.coefficients.lambda[i,q]=0; cmp.models.coefficients.Mn[i,q]=1
            warning(paste0('Reaction ',i,'-',q,' does not fit an exponential decay -- it was coerced to a horizontal line at binding saturation.'))
          }
        }
        cmp.models.data=list('FP'=c(data.scaled[,i,]),'tt'=rep(t,times=2)); try(cmp.models[[paste(i,'all',sep = '_')]]<-nls(FP~(1-Mn)*exp(-kn1*tt)+Mn,data=cmp.models.data,start=list('kn1'=1e-3,'Mn'=0),control=list(minFactor=1e-20,maxiter=1e2,warnOnly=TRUE)))
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
        plot(NULL,NULL,main = 'Comparison of Competition Reactions',xlab = 'Time (min)',ylab = 'Relative Polarization',ylim = c(0,1),xlim = c(min(t)/60,max(t)/60),cex.main=2,cex.axis=2,cex.lab=2)
      }
      if (scale.data==F){
        plot(NULL,NULL,main = 'Comparison of Competition Reactions',xlab = 'Time (min)',ylab = 'Polarization (mP)',ylim = c(min(na.omit(c(data))),max(na.omit(c(data)))),xlim = c(min(t)/60,max(t)/60),cex.main=2,cex.axis=2,cex.lab=2)
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
      par(fig=c(0,1,0,0.5),mar=c(5,6,3,1),oma = c(0,0,0,0),new=TRUE)
      plot(variant.concentrations*1e-3,rowMeans(k.off,na.rm = TRUE)*1e3,main = expression('Decoy-Dependence of k'[off]^obs),xlab = paste0('[Decoy] (µM)'),ylab = expression('k'['off']^obs*' x10'^-3*' (s'^-1*')'),type='p',cex.main=3,cex.lab=2.5,cex.axis=2.5)
      arrows(x0 = variant.concentrations*1e-3, x1 = variant.concentrations*1e-3, y0 = rowMeans(k.off,na.rm = TRUE)*1e3 - apply(k.off,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3, y1 = rowMeans(k.off,na.rm = TRUE)*1e3 + apply(k.off,MARGIN = c(1),FUN = sd,na.rm = TRUE)*1e3,code = 3,col = 'black',lwd = 1,angle = 90,length = 0.1)

      # Generate competition models
      fun.data=list('x'=(rep(variant.concentrations,times=4))[is.na(c(k.off-kt))==F],'y'=na.omit(c(k.off-kt)))
      if (estimate.initials==T){
        start.N=1
        start.ktheta=(rowMeans(k.off-kt,na.rm = TRUE)[1]-rowMeans(k.off-kt,na.rm = TRUE)[2])/(variant.concentrations[1]*1e-9-variant.concentrations[2]*1e-9)
        start.kn1=rowMeans(k.off-kt,na.rm = TRUE)[1]-(variant.concentrations[1]*1e-9)*start.ktheta
        start.Bhalf=variant.concentrations[which.min(abs(start.kn1/2-rowMeans(k.off-kt,na.rm = TRUE)))]
      }
      fun.model.opt2=nls(y~k.n1*(x*1e-9)^N/((x*1e-9)^N+(Bhalf*1e-9)^N)+k.theta*(x*1e-9),data = fun.data,start = list('k.n1'=start.kn1,'Bhalf'=start.Bhalf,'N'=start.N,'k.theta'=start.ktheta),control=list(minFactor=1e-10,maxiter=1e2,warnOnly=TRUE))
      if (match.optimizers==T & exists('fun.model.opt2')==T){
        fun.data<-list('x'=(rep(variant.concentrations,times=4))[is.na(c(k.off-kt))==F],'y'=na.omit(c(k.off-kt)),'Bhalf'=coefficients(fun.model.opt2)[['Bhalf']],'N'=coefficients(fun.model.opt2)[['N']])
        try(fun.model.opt<-nls(y~k.n1*(x*1e-9)^N/((x*1e-9)^N+(Bhalf*1e-9)^N),data = fun.data,start = list('k.n1'=coefficients(fun.model.opt2)[['k.n1']]),control=list(minFactor=1e-10,maxiter=1e2,warnOnly=TRUE)))
      }
      if (match.optimizers==F){
        fun.data=list('x'=(rep(variant.concentrations,times=4))[is.na(c(k.off-kt))==F],'y'=na.omit(c(k.off-kt)))
        try(fun.model.opt<-nls(y~k.n1*(x*1e-9)^N/((x*1e-9)^N+(Bhalf*1e-9)^N),data = fun.data,start = list('k.n1'=coefficients(fun.model.opt2)[['k.n1']],'Bhalf'=coefficients(fun.model.opt2)[['Bhalf']],'N'=coefficients(fun.model.opt2)[['N']]),control=list(minFactor=1e-10,maxiter=1e2,warnOnly=TRUE)))
      }
      if (exists('fun.model.opt')==T){
        lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(fun.model.opt,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt*1e3,col='green',lwd=3,lty='dashed')
        if (show.constants==T){
          cc.temp = signif((summary(fun.model.opt)[['coefficients']][1,1])*1e3,2); text(max(variant.concentrations)*1e-3*0.7,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.5+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=2.4,col='green')
        }
      }
      if (exists('fun.model.opt2')==T){
        lines((min(variant.concentrations):max(variant.concentrations))*1e-3,predict(fun.model.opt2,newdata=list('x'=min(variant.concentrations):max(variant.concentrations)))*1e3+kt*1e3,col='purple',lwd=3,lty='dashed')
        #legend(legend.location,legend = c('Classic Competition Model','Displacement-Transfer Model'),col = c('green','purple'),fill = c('green','purple'),cex = 3)
        legend(legend.location,legend = c('Classic Competition Model','Direct Transfer Model'),col = c('green','purple'),fill = c('green','purple'),cex = 3)
        if (show.constants==T){
          cc.temp = signif((summary(fun.model.opt2)[['coefficients']][1,1])*1e3,2); text(max(variant.concentrations)*1e-3*0.7,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.4+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[-1],' = ',cc.temp,'x10'^-3,' s'^-1,sep = "")),adj = c(0,0.5),cex=2.4,col='purple')
          cc.temp = signif((summary(fun.model.opt2)[['coefficients']][4,1]),2); text(max(variant.concentrations)*1e-3*0.7,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.3+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels = substitute(paste('k'[theta],' = ',cc.temp,' M'^-1,'s'^-1,sep = "")),adj = c(0,0.5),cex=2.4,col='purple')
        }
      }

      # Report and compare models
      if (exists('fun.model.opt')==T){
        print(paste('Classic Model Summary:',sep = ''),quote = F); print(summary(fun.model.opt),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (exists('fun.model.opt2')==T){
        print(paste('Direct Transfer Model Summary:',sep = ''),quote = F); print(summary(fun.model.opt2),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
      }
      if (exists('fun.model.opt2')==T & exists('fun.model.opt')==T){
        model.comparison=BIC(fun.model.opt,fun.model.opt2)
        print(paste("Models' BIC Statistics:",sep = ''),quote = F); print(model.comparison,quote = F)
        delta.BIC=BIC(fun.model.opt2)-BIC(fun.model.opt)
        print(paste('ΔBIC = ',delta.BIC,sep = ''),quote = F)
        print(paste(rep('#',times=100),collapse = ''),quote = F)
        if (delta.BIC<=-10){
          print('These data favor the DIRECT TRANSFER model!',quote = F)
          print(paste0('k.n1 = ',signif(summary(fun.model.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(fun.model.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
          print(paste0('k.theta = ',signif(summary(fun.model.opt2)[['coefficients']][4,1],2),' ± ',signif(summary(fun.model.opt2)[['coefficients']][4,2],2),' 1/M/s'),quote=F)
          print(paste(rep('#',times=100),collapse = ''),quote = F)
          if (show.constants==T){
            text(max(variant.concentrations)*1e-3*0.65,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.05+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels='BIC Favors >',cex = 5,col = 'purple',adj=c(1.5,0.5))
          }
        }
        if (delta.BIC>=5){
          print('These data favor the CLASSIC model!',quote = F)
          print(paste0('k.n1 = ',signif(summary(fun.model.opt)[['coefficients']][1,1],2),' ± ',signif(summary(fun.model.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
          print(paste(rep('#',times=100),collapse = ''),quote = F)
          if (show.constants==T){
            text(max(variant.concentrations)*1e-3*0.65,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.15+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels='BIC Favors >',cex = 5,col = 'green',adj=c(1.5,0.5))
          }
        }
        if (delta.BIC>-10 & delta.BIC<5){
          print('EITHER model may have produced these data!',quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
          print(paste0('k.n1 = ',signif(summary(fun.model.opt)[['coefficients']][1,1],2),' ± ',signif(summary(fun.model.opt)[['coefficients']][1,2],2),' 1/s'),quote=F)
          print('OR',quote=F)
          print(paste0('k.n1 = ',signif(summary(fun.model.opt2)[['coefficients']][1,1],2),' ± ',signif(summary(fun.model.opt2)[['coefficients']][1,2],2),' 1/s'),quote=F)
          print(paste0('k.theta = ',signif(summary(fun.model.opt2)[['coefficients']][4,1],2),' ± ',signif(summary(fun.model.opt2)[['coefficients']][4,2],2),' 1/M/s'),quote=F)
          print(paste(rep('#',times=100),collapse = ''),quote = F)
          if (show.constants==T){
            text(max(variant.concentrations)*1e-3*0.5,(diff(range(rowMeans(k.off,na.rm = TRUE)*1e3))*0.1+min(rowMeans(k.off,na.rm = TRUE)*1e3)),labels='BIC',cex = 5,col = 'grey',adj=c(1.5,0.5))
          }
        }
      }

      # Outlier Identification
      if (outliers[1]=='none' & default.mP.values==F){
        temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(fun.model.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]=='none' & default.mP.values==T){
        temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(fun.model.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=(((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==F){
        temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(fun.model.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))])[((colnames(raw.par)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
      if (outliers[1]!='none' & default.mP.values==T){
        temp.1=k.off; temp.1[is.na(temp.1)==FALSE]=(environment(fun.model.opt2[["m"]][["resid"]])[["resid"]])[is.na(temp.1)==FALSE]; outlyers=((((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))])[-((colnames(raw)[2:25])[c(seq(1,24,2),seq(2,24,2))] %in% outliers)==FALSE])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Obs.Num])[(EnvStats::rosnerTest(c(temp.1),k=5)$all.stats)$Outlier]
        if (length(outlyers)>=1){
          print(paste0('Outliers:'),quote = F); print(paste0(outlyers),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
        if (length(outlyers)==0){
          print(paste0('Outliers:'),quote = F); print(paste0('NONE'),quote = F); print(paste(rep('#',times=100),collapse = ''),quote = F)
        }
      }
    }
  }
  if (experiment.type=='STOICH'){
    nlm.Kd=coefficients(model.quadS)[1]
    nlm.S=coefficients(model.quadS)[2]
    nlm.Mn=coefficients(model.quadS)[3]
    nlm.Mx=coefficients(model.quadS)[4]
    show(summary(model.quadS))

    scale.dat.shi=(data.scaled.eq-nlm.Mn)/(nlm.Mx-nlm.Mn)
    nlm.sims=data.frame('y'=predict(model.quadS,newdata = list('E'=min(variant.concentrations):max(variant.concentrations))),'x'=min(variant.concentrations):max(variant.concentrations))
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

file.parse <- function(data.rows,background.rows,exp.type='full',par.file='par.csv',perp.file='perp.csv',path.to.file='./'){

  setwd(path.to.file)
  raw.par=read.csv(file = par.file,sep = ',',header = TRUE)
  raw.perp=read.csv(file = perp.file,sep = ',',header = TRUE)
  par.set=list(NULL)
  perp.set=list(NULL)
  if (exp.type=='single'){
    exp.n=length(data.rows)*2
    for (i in 1:exp.n){
      back.set=paste0(background.rows[i],1:24)
      if (i%%2==1){
        data.set=paste0(data.rows[round(i/2+0.01)],1:12)
      }
      if (i%%2==0){
        data.set=paste0(data.rows[round(i/2+0.01)],13:24)
      }
      whole.set=c('Time..s.',data.set,back.set)
      par.set[[i]]=raw.par[,(colnames(raw.par)%in%whole.set)]
      perp.set[[i]]=raw.perp[,(colnames(raw.perp)%in%whole.set)]
      write.table(x = par.set[[i]],file = paste0('par',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
      write.table(x = perp.set[[i]],file = paste0('perp',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
    }
  }
  if (exp.type=='half'){
    exp.n=length(data.rows)
    for (i in 1:exp.n){
      back.set=paste0(background.rows[i],1:24)
      data.set=paste0(data.rows[i],1:24)
      whole.set=c('Time..s.',data.set,back.set)
      par.set[[i]]=raw.par[,(colnames(raw.par)%in%whole.set)]
      perp.set[[i]]=raw.perp[,(colnames(raw.perp)%in%whole.set)]
      write.table(x = par.set[[i]],file = paste0('par',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
      write.table(x = perp.set[[i]],file = paste0('perp',i,'.txt'),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
    }
  }
  if (exp.type=='full'){
    exp.n=length(data.rows)/2
    for (i in 1:exp.n){
      back.set=paste0(background.rows[i],1:24)
      data.set=c(paste0(data.rows[2*(i-1)+1],1:24),paste0(data.rows[2*(i-1)+2],1:24))
      whole.set=c('Time..s.',data.set,back.set)
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

