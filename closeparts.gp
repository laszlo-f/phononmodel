reset
set border linewidth .5
set terminal postscript eps enhanced color dashed linewidth 3 size 2.5,8.*4/7 font "LMRoman12-Regular,14" fontfile "/usr/share/texmf/fonts/type1/public/lm/lmr12.pfb"
set pm3d map #interpolate 10,10
set key off
#unset surface
#set datafile missing "NaN"
hc=1239.84193  #eV nm
#set style data lines
#note: dgrid3d incompatable with NaN, apparently

set ylabel "Wavelength (nm)" 
set cblabel "Optical Density" rotate by -90
set zrange [-.01:0.02]
set yrange [589.38951:641.65454]
set xrange [-.3:4]
set tics out
set border lc rgb "#666666"
set tics textcolor rgb "black"

hc=1239.84193  #Si

hc2=0.000123984193 #(eV*cm)

k=8.6173324*10**(-5) #boltzmann eV/K

#Itoh et al 1975
g(x)=1/tanh(hc2*110/(2*k*x))
s(x)=hc2*(-207+5.2)*(g(x)-1)
temp=3.2
chirp(x,y)=x+1.528177e+03*(1/y-1/580.)


dynamics(x,start,tau,sigma)= (exp(sigma**2/(2*tau**2)-(x-start)/tau)*sqrt(pi/2)*sigma*(1-erf((sigma**2-tau*(x-start))/(sqrt(2)*sigma*tau))))
cadelay =    -7.562302e-02
caduration =  8.034848e-02
caamplitude = 1.836816e-03
sineduration= 4.994786e-02

ca(x)=caamplitude*sin((x-cadelay)/sineduration)*exp(-(x-cadelay)**2/(2*caduration**2))

wi(x,y)=6.650499e-02*dynamics(chirp(x,y),2.724239e-01,7.360879e-01,5.541777e-02)

mod=1.200434e-01
p1(x,y)=(((hc/y)-(2.033-.079+s(temp)))>0?mod*4.894297e-01*sqrt((hc/y)-(2.033-.079+s(temp)))*dynamics(chirp(x,y),2.724239e-01,3.073648e+01,5.541777e-02):0)
p2(x,y)=(((hc/y)-(2.033-.0138+s(temp)))>0?mod*sqrt((hc/y)-(2.033-.0138+s(temp)))*dynamics(chirp(x,y),2.724239e-01,9.155246e-01,5.541777e-02):0)

mod=1.200434e-01
ufod=4.835999e-01
enew=5.167266e-02
p3(x,y)=(((hc/y)-(2.033+enew+s(temp)))>0?ufod*mod*sqrt((hc/y)-(2.033+enew+s(temp)))*dynamics(chirp(x,y),2.724239e-01,9.818529e+00,5.541777e-02):0)

set output "partsc.eps"
set multiplot layout 4,1 offset .10,0.01
set colorbox user size .03,.15
TOP=.95
DY=.2
xsize=.70
ysize=.14*7/4
set ytics 20
set title offset -18,-2

#SET SAMPLES TO 500 FOR GOOD FILE
set samples 500
set isosamples 500
set title "(a)"# Coherent Artifact
set size xsize,ysize
set cbtics 0.001
set tmargin at screen TOP-0*DY
set bmargin at screen TOP-1*DY
set colorbox user origin .73,TOP-1*DY+.02
set xtics 1 out 
set xtics format ""
set ylabel offset 1
splot ca(chirp(x,y))

set title "(b)"# Wavelength-Independent
set size xsize,ysize
set tmargin at screen TOP-1*DY
set bmargin at screen TOP-2*DY
set colorbox user origin .73,TOP-2*DY+.02
set cbtics 0.002
unset xtics
set ylabel offset 0
splot wi(x,y)

set title "(c)" # Phonon A
set size xsize,ysize
set cbtics 0.001
set tmargin at screen TOP-2*DY
set bmargin at screen TOP-3*DY
set colorbox user origin .73,TOP-3*DY+.02
splot p2(x,y)

set title "(d)"# Phonon B
set ylabel offset 1
set size xsize,ysize
set cbtics 0.001
set xtics 1 out nomirror format "%g"
set tics textcolor rgb "black"
set xlabel "Probe Delay (ps)"
set tmargin at screen TOP-3*DY
set bmargin at screen TOP-4*DY
set colorbox user origin .73,TOP-4*DY+.02
splot p1(x,y)


unset multiplot

