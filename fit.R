library("scatterplot3d");
library(rgl);
library(plotrix);
options(digits=10);

hc=1239.84193  #eV*nm

hc2=0.000123984193 #(eV*cm)

k=8.6173324*10**(-5) #boltzmann eV/K

#used to compute exciton energy at temperature x
#Itoh et al 1975
g<-function(x){
	return(1/tanh(hc2*110/(2*k*x)));
}
s<-function(x){
	return(hc2*(-207+5.2)*(g(x)-1));
}

res<-500; #resolution
temp<-seq(0,1,length=res);
l<-seq(575,640,length=res);
t<-seq(-1,4,length=res);


# dynamics: 
# DSolve[{y'[x] == -y[x]/tau+pump*Exp[-x^2/(2*sigma^2)],y[a]==0}, y[x], x] 
# a is -Infinity

erf<-function(x){
	return(2*pnorm(x)-1);
}

dynamics<-function(x,start,tau,sigma){
	#pump=1
	return(exp(sigma^2/(2*tau^2)-(x-start)/tau)*sqrt(pi/2)*sigma*(1-erf((sigma^2-tau*(x-start))/(sqrt(2)*sigma*tau))));
}


#arguments
#
#lhold: wavelength
#th: time
#start: dynamics start time
#time1: lifetime of 79 meV phonon
#time2: lifetime of 13.8 meV phonon
#time3: lifetime of unkown peak
#time4: lifetime of wavelength independent decay
#timer: rise time
#cadelay: start time of coherent artifact
#caduration: how long artifact lasts
#sineduration: oscillation rate of coherent artifact; should be close to caduration
#mod: overall OD of absorption scale
#fod: OD of wavelength independent decay
#ufod: OD of unknown, high energy bump
#weakOD: relative OD of 79 meV phonon; should be about 0.165
#caamplitude: coherent artifact amplitude
#chirp: picosecond nanometers
#enew: threshold energy of unkown peak relative to 2.033 eV
af<-Vectorize(function(lhold,th,start,time1,time2,time3,time4,timer,cadelay,caduration,sineduration,mod,fod,ufod,weakOD,caamplitude,chirp,enew){
	temp<-3.2;#Kelvin
	thold<-  th+chirp*(1/lhold-1/580); #chirp
	value <-0;

	#absorption with two phonons + phonon created 
	if(((hc/lhold)-(2.033+enew+s(temp)))>0){
	       value<-
		       (ufod*sqrt((hc/lhold)-(2.033+enew+s(temp)))*dynamics(thold,start,time3,timer)#???
		       +sqrt((hc/lhold)-(2.033-.0138+s(temp)))*dynamics(thold,start,time2,timer)
	       +sqrt((hc/lhold)-(2.033-.079+s(temp)))*weakOD*dynamics(thold,start,time1,timer)
	       ) 
	#absorption with two phonons
	} else if(((hc/lhold)-(2.033-.0138+s(temp)))>0){
	       value<-
		       (sqrt((hc/lhold)-(2.033-.0138+s(temp)))*dynamics(thold,start,time2,timer)
	       +sqrt((hc/lhold)-(2.033-.079+s(temp)))*weakOD*dynamics(thold,start,time1,timer)#0.165
	       ) 

	#one phonon
	} else if(((hc/lhold)-(2.033-.079+s(temp)))>0){
	       value<-
		(
	       sqrt((hc/lhold)-(2.033-.079+s(temp)))*weakOD#0.165
	       ) *dynamics(thold,start,time1,timer);
	#no phonons
	} else {
		value <-0;
	}

	#wavelength independent background (???)
	value <-mod*value+fod*dynamics(thold,start,time4,timer);

	#coherentartifact
	#cadelay=.1;
	#caduration=.05;
	#sineduration#=.05;
	#caamplitude=0.005;#OD
	value<-value+caamplitude*sin((thold-cadelay)/sineduration)*exp(-(thold-cadelay)^2/(2*caduration^2));

	#convert absorption to OD
	return(value);
});

#make plot
#a<-outer(l,t,FUN="af");
#persp(l,t,z=a,ticktype="detailed",theta=120,phi=60,shade=.5,col="light green");

#read data
xyz<-read.table("allsub.xyz");
xyzroi<-subset(xyz,wavelength<(640)&time>-1.0);

#regression
print(date());
set.seed(1);
model=nls(ODsub ~ af(
		wavelength,time,start,time1,time2,time3,time4,timer,cadelay,caduration,sineduration,mod,fod,ufod,weakOD,caamplitude,chirp,enew
	),
	  data=xyzroi,
	  #data=xyzroi[sample(nrow(xyzroi),5e3),],
	  start=list(
		start=.254,time1=37,time2=0.9,time3=9,time4=.7,timer=0.0442,cadelay=-0.075,caduration=0.072,sineduration=0.049,mod=.194,fod=.057,ufod=0.477,weakOD=0.93,caamplitude=0.003,chirp=1305,enew=0.05678
	),control=c(tol=.00001)
	  ,trace=T);
print(date());
print(summary(model));

llong<-rep(l,100);
tlong<-sort(rep(t,100));
answer<-coefficients(model);
theory<-cbind(tlong,llong,af(llong,tlong,
			     answer[[1]],answer[[2]],answer[[3]],answer[[4]],answer[[5]],answer[[6]],answer[[7]],answer[[8]],answer[[9]],answer[[10]],answer[[11]],answer[[12]],answer[[13]],answer[[14]],answer[[15]],answer[[16]]));
colnames(theory)<-c("time","wavelength","Theory OD");
#open3d();
#plot3d(theory,size=1,col=color.scale(theory[,3],cs1=c(0,0,1,0,0),cs2=c(0,1,0,1,0),cs3=c(0,0,0,0,1)));

#open3d();
#plot3d(xyzroi,size=1,col=color.scale(xyzroi$ODsub,cs1=c(0,0,1,0,0),cs2=c(0,1,0,1,0),cs3=c(0,0,0,0,1)));
library("nlstools")
nlsContourRSS(model)
