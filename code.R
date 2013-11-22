## load of required libraries

library(raster)
library(solaR)

## introduce the working directory

dir<-as.character('/Users/usuario/Dropbox/ProyectoJavier/Programacion/Totales')
setwd(dir)

## definition of the daily temporal series. In this case 2005 daily time series (dTs)
## is defined with fBTd (solaR function)

dTs<-fBTd(mode='serie',start='01-01-2005',end='31-12-2005',format='%d-%m-%Y')

## definition of the intradaily temporal series (iTs) with as.POSIXct

iTs<-seq(as.POSIXct('2005-01-01 00:30:00'),as.POSIXct('2005-12-31 23:30:00'),by='hour')

## definition of month, day,hour, julian day and julian hour

dom<-c(31,28,31,30,31,30,31,31,30,31,30,31)
month<-rep(seq(1:12),dom*24)
df<-lapply(dom,function(x){t<-as.numeric(seq(1:x))})
day<-rep(do.call(c,df),each=24)
hour<-rep(1:24,365)
julianday<-rep(seq(1:365),each=24)
julianhour<-seq(1:8760)

## Solar field specifications

## reflectivity
LS3reflec<-0.94
## transmissivity
LS3trans<-0.955
## interceptation factor
LS3intf<-0.997
## absortivity
LS3abs<-0.955
## peak efficiency
nopico<-LS3reflec*LS3trans*LS3intf*LS3abs
## LS3-width
LS3width<-5.76
## LS3-focal distance
LS3focdis<-1.71
## LS3-SCA length
LS3length<-99
## LS3- operative surface
LS3oparea<-545
## Fouling factor
mirrorclean<-0.98
## HTF inlet temperature
HTFin<-293
## HTF outlet temperature
HTFout<-393
## Natural gas PCI [kJ/kg]
PCI<-39.900

## definition of stations where Totalccgt.csv is a file with the coordinates of stations.
ccgt<-read.csv2('Totalccgt.csv',sep=';',dec=',',header=TRUE)
est<-as.character(ccgt$Name[which(ccgt$Room==1)])

## Loop between 1 and 20 loops

for(n_loops in 1:20){
  modProd<-lapply(est,function(x){
    ## introduction of coordinates of the CCGT studied
    
    Lat<-ccgt$Latitude[which(ccgt$Room==1)][which(est==x)]
    Lon<-ccgt$Longitude[which(ccgt$Room==1)][which(est==x)]

    ## introduction of meteo data
    
    meteo<-read.csv2(paste(x,'Total.csv',sep=''))
    
    ## Apparent movement of the Sun from the Earth: incidence angle
    
    sol<-calcSol(Lat,local2Solar(iTs,Lon),sample='hour',EoT=TRUE,method='michalsky')
    incang<-r2d(acos(as.numeric(fTheta(sol,modeTrk='horiz')$cosTheta)[25:8784]))
    incang[which(is.na(incang))]<-0
    
    ## Modificator of incidence angle 
    
    modincang<-1-2.23073e-4*(incang)-1.1e-4*(incang^2)+3.18596e-6*(incang^3)-4.85509e-8*(incang^4)
    modincang[which(incang>80)]<-0
    
    ## Optical efficiency
    
    optef<-modincang*nopico
    
    ## Loss area per collector
    
    Lossarea<-LS3width*(LS3focdis+((LS3focdis*(LS3width^2))/(48*(LS3focdis^2))))*tan(d2r(incang))
    
    ## Heat loss from collector to environment (kWth per collector)
    
    Pcolenv<-(0.00154*(HTFm-meteo$TempMed)^2+0.2021*(HTFm-meteo$TempMed)-24.899+
                ((0.00036*(HTFm-meteo$TempMed)^2+0.2029*(HTFm-meteo$TempMed)+24.899)*(meteo$dni/900)*
                   cos(d2r(incang))))*LS3length/1000
    
    ## Potential thermal Power (kWth per collector)
    
    Psuncol<-LS3oparea*meteo$dni*cos(d2r(incang))/1000
    
    ## Thermal power from collector to fluid
    
    Pcolfluid<-(LS3oparea-Lossarea)*meteo$dni*cos(d2r(incang))*optef*mirrorclean/1000-Pcolenv
    Pcolfluid[which(Pcolfluid<0)]<-0
    Psolar_th<-n_loops*4*Pcolfluid ## 4 Solar Collection Assemblies per loop
    ef_st<-0.3425
    ef_parat<-0.94
    ef_process<-0.93
    Psolar_el<-Psolar_th*ef_st*ef_parat*ef_process/1000 ## in MW
    
    ## remove NA
    
    Psolar_el[which(is.na(Psolar_el))]<-0
    
    ## Annual thermal power from collector to fluid kWhth
    
    Pcolannual<-sum(Pcolfluid,na.rm=TRUE)
    
    ## Losses of power of combined cycle, gas turbine and steam turbine
    ## Supposing cummulative efficiency curves: GTrTaRHP (power)
    
    GTrTaRHP<-function(Ta,RH,P,Pb){
      efTa<-(-0.5024*Ta+107.536)/100
      efRH<-((1.05/90)*RH+99.3)/100
      efP<-((27/25)*((P/Pb)*100-100)+100)/100
      ef<-efTa*efRH*efP
      ef}
    
    ## Supposing cummulative efficiency curves: STrTaRHP (power)
    
    STrTaRHP<-function(Ta){
      efTa<-(6e-4*Ta^2-0.1579*Ta+102.24)/100
      ef<-efTa
      ef}
    
    ## Definition of inputs
    
    Ta<-meteo$TempMed
    RH<-meteo$HumMed
    P<-meteo$P
    DNI<-meteo$dni
    
    ## base pressure
    
    Pb<-((-27/2400)*ccgt$Altitude[which(ccgt$Name==x)]+100)/100*1013
    eGTrTa<-function(x){ef<-(-0.002*x^2-0.1237*x+102.28)}
    NG_cons_rate<-(20/260)*ccgt$TG[which(ccgt$Name==x)]
    efGTrTa<-eGTrTa(Ta)
    
    ### SCENARIO 1: solar boosting mode
    
    ## a) Natural gas consumption scenario 1
    
    NG_cons_sc1<-NG_cons_rate*efGTrTa/100
    
    ## b) GT Power output scenario 1
    
    Pgas_sc1<-GTrTaRHP(Ta,RH,P,Pb)*ccgt$TG[which(ccgt$Name==x)]
    
    ## c) Integrable solar power scenario 1
    
    Psol_int_sc1<-(-STrTaRHP(Ta)+1)*ccgt$TV[which(ccgt$Name==x)]
    
    ## remove negative values in Psol_int_sc1
    
    Psol_int_sc1[which(Psol_int_sc1<0)]<-0
    
    ## from integrable to energy integrated
    
    Psol_int_sc1i<-lapply(c(1:8760),function(t){
      if(Psolar_el[t]>Psol_int_sc1[t]){Psol_int_sc1[t]<-Psol_int_sc1[t]}else{Psol_int_sc1[t]<-Psolar_el[t]}
      Total<-Psol_int_sc1[t]})
    Psol_int_sc1<-do.call(c,Psol_int_sc1i)
    
    ## e) Dumping scenario 1
    
    Dumping_sc1<-Psolar_el-Psol_int_sc1
    Dumping_sc1[Dumping_sc1<0]<-0
    
    ## d) ST Power output scenario 1
    
    Pst_sc1<-STrTaRHP(Ta)*ccgt$TV[which(ccgt$Name==x)]+Psol_int_sc1
    
    ## d) Total CCGT power output scenario 1
    
    Pccgt_sc1<-Psol_int_sc1+Pst_sc1+Pgas_sc1
    
    ## f) Overall efficiency
    
    ef_ccgt_sc1<-Pccgt_sc1/(NG_cons_sc1*PCI)
    
    ### SCENARIO 2: solar dispatching mode
    
    ## a) Natural gas consumption scenario 2: calculation of reduced mass flow
    
    ratio_red<-(ccgt$TV[which(ccgt$Name==x)]-Psolar_el)/ccgt$TV[which(ccgt$Name==x)]*100
    ratio_red[which(is.na(ratio_red))]<-100
    
    ## relationship between load of GT and efficiency
    
    eGTrLoad<-function(load){ef<-(-0.0058*load^2+1.3125*load+26.645)/100} 
    efGTrLoad<-eGTrLoad(ratio_red)
    efGTrLoad[efGTrLoad==0.9989500]<-1
    NG_cons_sc2<-NG_cons_rate*efGTrTa*efGTrLoad/100
    
    ## b) GT power output scenario 2
    
    Pgas_sc2<-efGTrLoad*GTrTaRHP(Ta,RH,P,Pb)*ccgt$TG[which(ccgt$Name==x)]
    
    ## c) ST power output scenario 2
    
    Pst_sc2<-rep(ccgt$TV[which(ccgt$Name==x)],8760)
    Pst_sc2i<-lapply(c(1:8760),function(t){
      if(Psolar_el[t]>0){Pst_sc2[t]<-ccgt$TV[which(ccgt$Name==x)]}else{Pst_sc2[t]<-ccgt$TV[which(ccgt$Name==x)]*STrTaRHP(Ta[t])}
      Total<-Pst_sc2[t]})
    Pst_sc2<-do.call(c,Pst_sc2i)
    
    ## d) CCGT power output scenario 2
    
    Pccgt_sc2<-Pgas_sc2+Pst_sc2
    
    ## e) Dumping
    
    Dumping_sc2<-rep(0,8760)
    
    ## f) Solar power integration
    
    Psol_int_sc2<-Psolar_el
    
    ## g) Overall efficiency
    
    ef_ccgt_sc2<-Pccgt_sc2/(NG_cons_sc2*PCI)
    
    ### SCENARIO 0: CCGT
    
    ## a) Natural gas consumption
    
    NG_cons_sc0<-NG_cons_rate*efGTrTa/100
    
    ## b) GT power output scenario 0
    
    Pgas_sc0<-GTrTaRHP(Ta,RH,P,Pb)*ccgt$TG[which(ccgt$Name==x)]
    
    ## c) ST power output sceneario 0
    
    Pst_sc0<-STrTaRHP(Ta)*ccgt$TV[which(ccgt$Name==x)]
    
    ## d) CCGT power output scenario 0
    
    Pccgt_sc0<-Pgas_sc0+Pst_sc0
    
    ## e) Solar power integration
    
    Psol_int_sc0<-rep(0,8760)
    
    ## g) Overall efficiency
    
    ef_ccgt_sc0<-Pccgt_sc0/(NG_cons_sc0*PCI)
    
    Total<-data.frame(iTs,month,day,julianday,hour,julianhour,DNI,Ta,RH,P,
                      incang,modincang,optef,Lossarea,Pcolenv,Psuncol,Pcolfluid,Psolar_th,Psolar_el,efGTrTa,
                      NG_cons_sc1,Pgas_sc1,Psol_int_sc1,Dumping_sc1,Pccgt_sc1,Pst_sc1,ef_ccgt_sc1,ratio_red,efGTrLoad,
                      NG_cons_sc2,Pgas_sc2,Pst_sc2,Pccgt_sc2,Dumping_sc2,Psol_int_sc2,ef_ccgt_sc2,NG_cons_sc0,
                      Pgas_sc0,Pst_sc0,Pccgt_sc0,Psol_int_sc0,ef_ccgt_sc0)
    
    Total})
  old<-getwd()
  save(modProd,file=paste('Todo',n_loops,'.RData',sep=''))
  setwd(old)
}

## Annual results for scenarios 0, 1 and 2

Annual<-lapply(c(1:20),function(x){
  load(paste('Todo',x,'.RData',sep=''))
  t<- lapply(modProd,function(x){
    Pccgt_sc0<-sum(x$Pccgt_sc0)
    Pccgt_sc1<-sum(x$Pgas_sc1)+sum(x$Pst_sc1)+sum(x$Psol_int_sc1)
    Pccgt_sc2<-sum(x$Pccgt_sc2)
    Psol_int_sc1<-sum(x$Psol_int_sc1)
    Psol_int_sc2<-sum(x$Psol_int_sc2)
    Dumping_sc1<-sum(x$Dumping_sc1)
    NG_cons_sc0<-sum(x$NG_cons_sc0)
    NG_cons_sc2<-sum(x$NG_cons_sc2)
    Pgas_sc1<-sum(x$Pgas_sc1)
    Pst_sc1<-sum(x$Pst_sc1)
    Pgas_sc2<-sum(x$Pgas_sc2)
    Pst_sc2<-sum(x$Pst_sc2)
    Pst_sc0<-sum(x$Pst_sc0)
    Pgt_sc0<-sum(x$Pgas_sc0)
    Total<-c(Pccgt_sc0,Pst_sc0,Pgt_sc0,Pccgt_sc1,Pccgt_sc2,Psol_int_sc1,Psol_int_sc2,Dumping_sc1,NG_cons_sc0,NG_cons_sc2,
             Pgas_sc1,Pgas_sc2,Pst_sc1,Pst_sc2)
  })
  Total<-data.frame(do.call(rbind,t))
  names(Total)<-c('Pccgt_sc0','Pst_sc0','Pgt_sc0','Pccgt_sc1','Pccgt_sc2','Psol_int_sc1','Psol_int_sc2','Dumping_sc1','NG_cons_sc0',
                  'NG_cons_sc2','Pgas_sc1','Pgas_sc2','Pst_sc1','Pst_sc2')
  return(Total)
})
save(Annual,file='Annual.RData')

## Annual results for scenarios 3, 4 and 5

horas_operacion<-seq(11,17,1)
horas<-lapply(horas_operacion,function(x){
  t<-which(modProd[[1]]$hour==x)
})
horas<-sort(do.call(c,horas))

Annual_sc34<-lapply(c(1:20),function(x){
  load(paste('Todo',x,'.RData',sep=''))
  t<- lapply(modProd,function(x){
    Pccgt_sc0<-sum(x$Pccgt_sc0[horas])
    Pccgt_sc3<-sum(x$Pgas_sc1[horas])+sum(x$Pst_sc1[horas])+sum(x$Psol_int_sc1[horas])
    Pccgt_sc4<-sum(x$Pccgt_sc2[horas])
    Psol_int_sc3<-sum(x$Psol_int_sc1[horas])
    Psol_int_sc4<-sum(x$Psol_int_sc2[horas])
    Dumping_sc3<-sum(x$Dumping_sc1[horas])
    NG_cons_sc0<-sum(x$NG_cons_sc0[horas])
    NG_cons_sc4<-sum(x$NG_cons_sc2[horas])
    Pgas_sc3<-sum(x$Pgas_sc1[horas])
    Pst_sc3<-sum(x$Pst_sc1[horas])
    Pgas_sc4<-sum(x$Pgas_sc2[horas])
    Pst_sc4<-sum(x$Pst_sc2[horas])
    Pst_sc0<-sum(x$Pst_sc0[horas])
    Pgt_sc0<-sum(x$Pgas_sc0[horas])
    Total<-c(Pccgt_sc0,Pst_sc0,Pgt_sc0,Pccgt_sc3,Pccgt_sc4,Psol_int_sc3,Psol_int_sc4,Dumping_sc3,NG_cons_sc0,NG_cons_sc4,
             Pgas_sc3,Pgas_sc4,Pst_sc3,Pst_sc4)
  })
  Total<-data.frame(do.call(rbind,t))
  names(Total)<-c('Pccgt_sc5','Pst_sc5','Pgt_sc5','Pccgt_sc3','Pccgt_sc4','Psol_int_sc3','Psol_int_sc4','Dumping_sc3','NG_cons_sc00',
                  'NG_cons_sc4','Pgas_sc3','Pgas_sc4','Pst_sc3','Pst_sc4')
  return(Total)
})
save(Annual_sc34,file='Annual_sc34.RData')
