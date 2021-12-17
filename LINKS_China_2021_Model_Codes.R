install.packages("deSolve") ##Can delete this code after installing "desolve"
library("deSolve")


## MODELING PART
model <- function (t, x, params) {
  RS <- x[1] ##Riboswitch
  RSWoff <- x[2] ##Off-state riboswitch
  RSWon <- x[3] ##On-state riboswitch
  cI <- x[4] ##cI monomers
  cI2 <- x[5] ##cI dimers
  mRNA <- x[6] ##mRNA of tnaA-FL-FMO
  output <- x[7] ##protein of tnaA-FL-FMO
  W <- x[8] ##trp concentration
  diW <- x[9] ##6-Br-trp concentration
  indigo <- x[10] ##Indigo concentration
  dixindigo <- x[11] ##Tyrian purple concentration
  
  
  DNA <- 12 ##Riboswitch plasmid copy number
  DNA2 <- 12 ##tnaA-FL-FMO plasmid copy number
  kTX <- params["kTX"] ##Transcription rate of constitutive promoter
  kDegm <- params["kDegm"] ##mRNA degradation rate 
  k1 <- params["k1"] ##Association rate between trp and riboswitch
  k1r <- params["k1r"] ##Dissociation rate between trp and riboswitch
  kon <- params["kon"] ##Rate of ribosiwtch turning on
  koff <- params["koff"] ##Rate of ribosiwtch turning off  
  kTL <- params["kTL"] ## Translation rate of cI
  kLk <- params["kLk"] ##Leakage translation rate on off-state riboswitch 
  kDegp <- params["kDegp"] ##Protein degradation rate
  k2 <- params["k2"]
  k2r <- params["k2r"]
  kTX2 <- params["kTX2"] ##Transcription rate of pRP when cI is not present at all
  kR <- params["kR"] ##Repression constant of cI repressor 
  n <- params["n"] ##Hill coefficient of cI 
  kTL2<-params["kTL2"] ##Translation rate of tnaA-FL-FMO
  Vmax <- params["Vmax"] ##Vmax of fre-sttH
  Km <- params["Km"]  ##Km of fre-sttH
  kb <- params["kb"] ##Binding constant of cI2 and DNA2
  kbr <- params["kbr"] ##Dissociation constant of cI2 and DNA2
  kcat <- params["kcat"] ##Rate constant of trp-tnaA-FL-FMO complex dissociating to give indigo

  
  #Equations
  dRSdt <- kTX*DNA - kDegm*RS + k1r*RSWoff
  dRSWoffdt <- k1*RS*W - k1r*RSWoff 
  dRSWondt <- kon*RSWoff - koff*RSWon
  dcIdt <- kTL*RSWon + kLk*(RSWoff+RS) - kDegp*cI
  dcI2dt <- k2*cI^2 - kDegp*cI2
  dmRNAdt <- DNA2*kTX2*1/(1+(cI2/kR)^n) - kDegm*mRNA
  doutputdt <- kTL2*mRNA - kDegp*output
  dWdt<- - k1*RS*W + k1r*RSWoff  - (Vmax*W)/(Km + W) - (kcat*W*output)/(0.01 + W)
  ddiWdt <- (Vmax*W)/(Km + W) - (kcat*diW*output)/(0.01 + diW)
  dindigodt <- (kcat*W*output)/(0.01 + W)
  ddixindigodt <- (kcat*diW*output)/(0.01 + diW)
  
  
  dxdt <- c(dRSdt, dRSWoffdt, dRSWondt, dcIdt, dcI2dt, dmRNAdt, doutputdt, dWdt, ddiWdt, dindigodt, ddixindigodt)
  list(dxdt)
}

  #Parameters
params <- c(kTX=0.2229, kDegm=5.5*10^-3, k1=2.16*10^-6, k1r=7.99*10^-3,
            kon=1*10^-1, koff=1*10^-5,  kTL=0.5, kLk=0.0005, kDegp=1.16*10^-3, k2=10^-6, k2r=0.0010,
            kTX2=0.05912, kTL2=0.0213, kR=4*10, n=2, Vmax=7.08778*10^-8, Km=4.060915*10^-3, kb=0.0006, kbr=0.1, kcat=10^-10.5)
xstart<- c(RS=0, RSWoff=0, RSWon=0, cI=0, cI2=0,mRNA=0, output=0, W=2.5*10^-3 , diW=0, indigo=0, dixindigo=0)

times<-seq(from=0, to=1004800, by=1000)
out<-ode(
  func = model,
  y = xstart,
  times = times,
  parms = params
)
out.df <- as.data.frame(out)









## GRAPH PLOTTING PART (ggplot codes)
library(tidyverse)
out.df<-out.df%>%
  mutate('TnaA-L-Fmo' =output)%>%
  mutate('Trp'=W)%>%
  mutate('6-Br-Trp'=diW)%>%
  mutate('Indigo'=indigo)%>%
  mutate('Tyrian purple'=dixindigo)%>%
  mutate('Ratio'=dixindigo/indigo)
library(ggplot2)

library(tidyr)
#1.Protein plot
df <- gather(out.df, key = Protein, value = Amount, 
             c("cI", "cI2", "TnaA-L-Fmo"))

ggplot(df, aes(x=time, y = Amount, group = Protein, colour = Protein)) + 
  geom_line()+
  xlab('Time(s)')+
  ylab('Amount of cI, cI2 and TnaA-L-Fmo')

##Note: remove each individual protein in gather function and 
##in y label to get plot for each individual protein

#2.Amino acid plot
df2 <- gather(out.df, key = measure, value = Concentration, 
              c("Trp", "6-Br-Trp"))

ggplot(df2, aes(x=time, y = Concentration, group = measure, colour = measure)) + 
  geom_line()+
  xlab('Time(s)')+
  ylab('Concentration(M)')+
  labs(color = "Amino acid")

#3.Indigo and Tyrian purple plot
df3 <- gather(out.df, key = measure, value = indigo, 
              c("Indigo", "Tyrian purple"))

ggplot(df3, aes(x=time, y = indigo, group = measure, colour = measure)) + 
  geom_line()+
  xlab('Time(s)')+
  ylab('Concentration(M)')+
  labs(color = "Dye")
#4. Ratio
ggplot(out.df, aes(x=time))+
  geom_line(aes(y=Ratio), col='blue')+
  xlab('Time(s)')+
  ylab('Tyrian purple : Indigo ratio')