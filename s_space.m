function y=s_space(xstart,xend,num)
% S_SPACE(XSTART,XEND,NUM) generates num evenly spaced sampled 
% values from xstart to (xend-dx), where dx is the sample 
% spacing. This is useful in FFT analysis where we generate 
% sampled periodic functions. Example: generate 1000
% sampled frequencies from 0 to 100MHz via f=s_space(0,100,1000); 
% In this case the last value of f will be 99.9 MHz and the 
% sampling interval will be 100/1000= 0.1 MHz.
% script written by Les Schmerr, Iowa Univ. 
% in google: s_space.m matlab script schmerr
% see info in : http://www.public.iastate.edu/~e_m.350/FFT%205.pdf
%
ye=linspace(xstart,xend,num+1);
y=ye(1:num);
