function [x,noi]=afnlm(a1,f,s,l1,l2,fmi,fmx,sigma,decon,deconnoi)
% afnlm: adptive frequency NLM for the noisy data
% 
% BY Hang Wang, April, 2021
% INPUT
% a1: input noisy data
% f:  the half side length of the searching area
% s:  the half side length of the patch area
% l1,l2: filtering parameter
% fmi,fmx: frequency bandwidth in sample points
% sigma: variance of gauss kernal
% decon: denoised data by the pre-proecssing
% deconnoi:  noise data by the pre-proecssing

% OUTPUT
% x:  denoised data
% noi: removed noise

[m,~]=size(a1);
nfa=fft(a1);
nfar=real(nfa);
nfai=imag(nfa);

deconfr=real(fft(decon));
deconfi=imag(fft(decon));
deconnoifr=real(fft(deconnoi));
deconnoifi=imag(fft(deconnoi));

tic;
ter=nfar(fmi:fmx,:);%real part
z=(ter-min(min(ter)))./(max(max(ter))-min(min(ter)))*255;
[xr,~]=nlm_f(z,deconfr,deconnoifr,f,s,sigma,l1,l2,fmi,fmx);
xr=xr.*(max(max(ter))-min(min(ter)))/255+min(min(ter));

tei=nfai(fmi:fmx,:);%imagery part
z=(tei-min(min(tei)))./(max(max(tei))-min(min(tei)))*255;
[xi,~]=nlm_f(z,deconfi,deconnoifi,f,s,sigma,l1,l2,fmi,fmx);
xi=xi.*(max(max(tei))-min(min(tei)))/255+min(min(tei));

Xr=a1*0;%Conjugate symmetry
Xr(1:fmx,:)=[nfar(1:fmi-1,:);xr];
for i=2:floor((m-2)/2+2)
    Xr(m-i+2,:)=Xr(i,:);
end

Xi=a1*0;%Conjugate symmetry
Xi(1:fmx,:)=[nfai(1:fmi-1,:);xi];
for i=2:floor((m-2)/2+2)
    Xi(m-i+2,:)=-Xi(i,:);
end

xr1=Xr;
xi1=Xi;
X=xr1+1i.*xi1;
x=real(ifft(X));
noi=a1-x;

toc;
end