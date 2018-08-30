

load("~/Science/Mammals Trophic Structure/EFWOGE/Environments/1-3MB.9TG.RData")

# solo AMD
load("~/Science/Mammals Trophic Structure/EFWOGE/Environments/1.MB.AMD_TG.RData")

library(e1071)

rm(list=ls(all=TRUE))

setwd("~/Science/Mammals Trophic Structure/Mambirds/Data")
MBData <- read.table(file="MBData.csv",header=T,sep=",")
MBData <- MBData[-c(1)] # Quita X

spc<-as.data.frame(MBData$Spp)
setwd("~/Science/Mammals Trophic Structure/Mambirds/Data")
# write.csv(spc, 'spc.csv') 

MBData <- MBData[-c(1)] # Quita Spp

# ---------------         Number of TG

its<-100 # n? de veces que se repiten los an?lisis, para evitar el efecto del azar sobre el clustering

nin<-2
nsp<-1000
r<-nsp-nin

nc<-ncol(MBData)
nr<-nrow(MBData)

Maxpms<-c()

for(k in 1:its){ # para que repita los an?lisis its veces (para evitar el efecto del azar sobre el clustering)
  
  Maxpm<-c()
  
  for(i in nin:nsp){ # para que repita los an?lisis para cada n? de clusters
    
    cl<-cmeans(MBData,i,20,verbose=FALSE,method="cmeans",m=2)  # hace fuzzy clust. con distinto n? de grupos, desde nin a nsp.
    
    Prob<-cl$membership # matriz con las probabilidades de pertenencia de cada muestra a cada grupo
    
    # Calcula la probabilidad de pertenencia de cada muestra al grupo al que fue asignada en "cluster".
    Maxprob<-c()
    for(j in 1:nr){
      Maxprob[[j]] <- c(max(Prob[j,]))
    }
    
    Mpm<-mean(Maxprob)-1/i # Calcula la media de las probabilidades max. calculadas en el paso anterior.
    Maxpm<-cbind(Maxpm,Mpm) # crea una BD con las Mpm de cada uno de los an?lisis, con distinto n? de grupos, de nin a nsp.
    N<-ncol(Maxpm)
    colnames(Maxpm)[N]<-paste(i,"Cls",sep="") # nombra el n? de clusters al que corresponde cada valor de CRc
    
  }
  
  Maxpms<-rbind(Maxpms,Maxpm) # crea un BD con los valores de Maxpm de cada itr
  
  if (i==nsp) print (k) # para que vaya indicando, mientras corre el script, por qu? n? de itr va
}

colMax <- function (colData) {apply(colData, MARGIN=c(2), max)}

Maxpmean<-data.frame(colMax(Maxpms)) # calcula el max. de todos los valores de CRs obtenidos para cada cada n? de clusters
colnames(Maxpmean)[1]<-"Maxpmean"

clusters<-c(nin:nsp)
Results<-cbind(clusters,Maxpmean)

windows()
plot(nin:nsp,Results$Maxpmean, xlab="No of Clst",ylab="Prob (Max) mean") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = (ls()[ls() != "MBData"]))

N<-9

its<-200 # n? de veces que se repiten los an?lisis
nr<-nrow(MBData)

Maxpm<-c(0)
Clst<-c()
MaxProbs<-c()
matprobs <- matrix(nrow=nr,ncol=N)


for(i in 1:its){ # para que repita los an?lisis its veces (para evitar el efecto del azar sobre el clustering)
  
  
  cl<-cmeans(MBData,N,20,verbose=FALSE,method="cmeans",m=2)  # hace fuzzy clust. con distinto n? de grupos, desde nin a nsp.
  
  cluster<-as.numeric(cl$cluster)# clasifica las muestras en los i grupos
  Prob<-cl$membership # matriz con las probabilidades de pertenencia de cada muestra a cada grupo
  
  Maxprob<-c()
  for(j in 1:nr){
    Maxprob[[j]] <- c(max(Prob[j,]))
  }
  
  Mpm<-mean(Maxprob)-1/N
  
  cat (i, "\n")
  
  if (Mpm>Maxpm){
    Maxpm<-Mpm
    Clst<-cluster
    MaxProbs<-Maxprob
    matprobs <- Prob
    print(Maxpm)
  }
}

matprobs <- round(matprobs,3)
MaxProbs <- round(MaxProbs,3)

# Clst es el resultado IMP del loop


Clst[Clst=="1"] <- "Nct"
Clst[Clst=="2"] <- "Grn"
Clst[Clst=="3"] <- "IFd"
Clst[Clst=="4"] <- "GIn"
Clst[Clst=="5"] <- "Pln"
Clst[Clst=="6"] <- "Frg"
Clst[Clst=="7"] <- "Crn"
Clst[Clst=="8"] <- "HIn"
Clst[Clst=="9"] <- "IFr"


MBDataC<-cbind(MBData,Clst)
MBDataC$Clst<-as.factor(MBDataC$Clst)
MBDataC <- MBDataC[-c(1)] # Quita Spp

setwd("~/Science/Mammals Trophic Structure/Mambirds/Data")
# write.csv(MBDataC, 'MBDataC.csv') 
Clst<-as.data.frame(Clst)
# write.csv(Clst, 'Clst.csv')


SppTG<-as.data.frame(cbind(spc,Clst))

colnames(SppTG)[1] <- "Spp"
#write.csv(SppTG, 'SppTG.csv') 


#------

Count<-table(MBDataC$Clst)
Count
TG<-unique(MBDataC$Clst)

windows()
par(mfrow=c(3,3), mar=c(4,3,2,2.2))

for(i in TG){
  s<-subset(MBDataC,subset=Clst==i) 
  boxplot(s[,1:10],outline=F, range=1,ylim=c(0,100))
  title(paste(i,"(",Count[i],")"), line = -1.5)
}
