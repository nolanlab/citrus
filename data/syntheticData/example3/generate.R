library("flowCore")

dataDir = "Desktop/work/citrus/data/syntheticData/example2/"


genBasePop = function(x){
  pop1 = matrix(rnorm(nEvents*percents[1]*2),ncol=2)
  pop2 = matrix(rnorm(nEvents*percents[2]*2,mean=3.5),ncol=2)
  pop3 = cbind(rnorm(nEvents*percents[3],mean=3.5),rnorm(nEvents*percents[3]))
  return(rbind(pop1,pop2,pop3))
}


nEvents = 10000
percents = c(0.6,.25,.15)
for (i in 1:10){
    for (diseaseState in c("healthy","diseased")){
      trueLabels = c(rep(1,nEvents*percents[1]),rep(2,nEvents*percents[2]),rep(3,nEvents*percents[3]))      
      fBase = rnorm(1,mean=1)
      if (diseaseState=="diseased"){
        fStim = fBase + abs(rnorm(1,sd=.5,mean=1))  
      } else {
        fStim = fBase
      }
      print(paste(fBase,fStim))
      
      f1Unstim = rnorm(nEvents,mean=2)
      f1Stim1 = rnorm(nEvents,mean=2)
      f2Unstim = rnorm(nEvents,mean=fBase)
      f2Stim1 = rnorm(nEvents,mean=fBase)
      f2Stim1[trueLabels==3]=rnorm(sum(trueLabels==3),mean=fStim)
    
      if (diseaseState=="diseased"){
        patientId = i+10;  
      } else {
        patientId=i;
      }
      
      
      patientIdStim = paste("Patient",sprintf("%02.0f",patientId),"_",diseaseState,"_stim1.fcs",sep="")
      patientIdUnstim = paste("Patient",sprintf("%02.0f",patientId),"_",diseaseState,"_unstim.fcs",sep="")
      patientUnstimData =  cbind(genBasePop(),f1Unstim,f2Unstim,trueLabels)
      patientStimData =  cbind(genBasePop(),f1Stim1,f2Stim1,trueLabels)
      patientUnstimData = patientUnstimData[sample(1:nEvents),]
      patientStimData = patientStimData[sample(1:nEvents),]
      
      colnames(patientUnstimData) = c("LineageMarker1","LineageMarker2","FunctionalMarker1","FunctionalMarker2","trueLabel")
      colnames(patientStimData) = c("LineageMarker1","LineageMarker2","FunctionalMarker1","FunctionalMarker2","trueLabel")
      
      plot(density(patientUnstimData[patientUnstimData[,"trueLabel"]==2,3]))
      lines(density(patientStimData[patientStimData[,"trueLabel"]==2,3]))
      
      plot(density(patientUnstimData[patientUnstimData[,"trueLabel"]==3,3]))
      lines(density(patientStimData[patientStimData[,"trueLabel"]==3,3]))
      
      plot(density(patientUnstimData[patientUnstimData[,"trueLabel"]==2,4]))
      lines(density(patientStimData[patientStimData[,"trueLabel"]==2,4]))
      
      plot(density(patientUnstimData[patientUnstimData[,"trueLabel"]==3,4]))
      lines(density(patientStimData[patientStimData[,"trueLabel"]==3,4]))      
      
      write.FCS(flowFrame(patientUnstimData),filename=file.path(dataDir,patientIdUnstim))
      write.FCS(flowFrame(patientStimData),filename=file.path(dataDir,patientIdStim))
  }
  
  
}
