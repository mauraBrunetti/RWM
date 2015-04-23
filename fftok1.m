function [a,phi]=fftok1(x,Ti,fa,df) 
t=(0:(length(x)-1))/fa;
xc=x.*exp(2*1i*pi.*t./Ti).';
xc=sum(xc(1:end-1)+xc(2:end))/fa*df;      		  
a=abs(xc); phi=angle(xc); 
