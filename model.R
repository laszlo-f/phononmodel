library("scatterplot3d");
library(rgl);
library(plotrix);
options(digits=10);

hc=1239.84193  #Si

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

dynamics<-function(x,tau,pump,sigma){
	return(exp(sigma^2/(2*tau^2)-x/tau)*sqrt(pi/2)*pump*sigma*(1-erf((sigma^2-tau*x)/(sqrt(2)*sigma*tau))));
}


af<-Vectorize(function(lhold,th){
	temp<-4;#Kelvin
	delay<-.4
	thold<-th-delay+1000*(1/lhold-1/580); #chirp
	value <-0;

	#absorption with two phonons + phonon created 
	if(((hc/lhold)-(2.033+.0138+s(temp)))>0){
	       value<-
		       (.5*sqrt((hc/lhold)-(2.033+.0138+s(temp)))#???
		       +sqrt((hc/lhold)-(2.033-.0138+s(temp)))
	       +sqrt((hc/lhold)-(2.033-.079+s(temp)))*0.165
	       ) *dynamics(thold,2,1,.05);
	#absorption with two phonons
	} else if(((hc/lhold)-(2.033-.0138+s(temp)))>0){
	       value<-
		       (sqrt((hc/lhold)-(2.033-.0138+s(temp)))
	       +sqrt((hc/lhold)-(2.033-.079+s(temp)))*0.165
	       ) *dynamics(thold,2,1,.05);

	#one phonon
	} else if(((hc/lhold)-(2.033-.079+s(temp)))>0){
	       value<-
		(
	       sqrt((hc/lhold)-(2.033-.079+s(temp)))*0.165
	       ) *dynamics(thold,2,1,.05);
	#no phonons
	} else {
		value <-0;
	}

	#wavelength independent background (???)
	value <-.2*value+.05*dynamics(thold,2,1,.05);

	#coherentartifact
	cadelay=.03;
	caduration=.025;
	sineduration=.025;
	caamplitude=0.005;#OD
	value<-value+caamplitude*sin((thold-cadelay)/sineduration)*exp(-(thold-cadelay)^2/(2*caduration^2));

	#convert absorption to OD
	return(value);
});

#make plot
#a<-outer(l,t,FUN="af");
#persp(l,t,z=a,ticktype="detailed",theta=120,phi=60,shade=.5,col="light green");

llong<-rep(l,100);
tlong<-sort(rep(t,100));
theory<-cbind(tlong,llong,af(llong,tlong));
colnames(theory)<-c("time","wavelength","Theory OD");
open3d();
plot3d(theory,size=1,col=color.scale(theory[,3],cs1=c(0,0,1,0,0),cs2=c(0,1,0,1,0),cs3=c(0,0,0,0,1)));

#read data
xyz<-read.table("allsub.xyz");
xyzroi<-subset(xyz,wavelength<(640)&time>-1.0);
open3d();
plot3d(xyzroi,size=1,col=color.scale(xyzroi$ODsub,cs1=c(0,0,1,0,0),cs2=c(0,1,0,1,0),cs3=c(0,0,0,0,1)));
