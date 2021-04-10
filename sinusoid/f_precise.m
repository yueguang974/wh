function [n,ffs]=f_precise(p,m,ninter,dt)
% precise frequency analysis
% 
% BY Hang Wang, April, 2021
% INPUT
% p: input noisy data
% m: length before the interpolation
% ninter:  the number of the interpolation points
% dt:  time interval


% OUTPUT
% n:  precise period (sample points)
% ffs: precise frequency

pp=abs(fft(p));
Fp=1/dt;
fmax=(find(pp(1:end/2,1)==max(pp(1:end/2,1)))-1)*Fp/m;
%rough dominant frequency

ffs=precise_dft(p,floor(fmax)-2,floor(fmax)+2,400,ninter,dt);
%search within the scope of [f-2,f+2], divide it into 400 small intervals

n=floor(Fp/ffs*ninter);
end