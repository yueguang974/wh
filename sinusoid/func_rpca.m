function [p4,ffs]= func_rpca(sn,ninter,lr,rk,wid,rto,rd,dt)
% randomized PCA for the sinusoidal noise attenuation
% 
% BY Hang Wang, April, 2021
% INPUT
% sn: input noisy data
% ninter:  the number of the interpolation points
% lr:  the number of periods used for the RPCA matrix (lateral size of C)
% rk:  truncted rank value of C
% wid: bandpass filter radius in sample-points unit
% rto: scaling factor for the ranks
% rd: if using randomization 0=no, 1=yes
% dt: time interval


% OUTPUT
% p4:  denoised data
% ffs: removed noise

n=length(sn);
p=interp1(1:n,sn,1:(1/ninter):n,'spline')';%interpolation of noisy data


m=size(p(1:ninter:end,:),1);%length before the interpolation
[n,ffs]=f_precise(p,m,ninter,dt);%calculate period after the interp
NN=size(p,1);%length after the interpolation
N=n*floor(NN/n);%how many complete periods do the data contain(not count the surplus data points at the end)


po=[p(1:N,1);p(N+1:NN,1);p(NN-n+1:N,1)];
%Add some data at the end to make its length the integral multiple of the period

p=po;
kp=cat(1,repmat(p(1:n,:),lr/2,1),p,repmat(p(end-n+1:end,:),lr/2,1));%
po=kp;
%padding at the begining and the end with lr/2 periods




NN=NN+lr*n;
NNN=size(po,1);
%length after the padding


cen=NNN/n;
pf=fft(po);
fil=ones(NNN,1);
if cen-wid<=1
    wid=cen-1;
end
fil(cen-wid:cen+wid)=0;

invf=zeros(NNN,1);
invf(1:floor((2+NNN)/2))=pf(1:floor((2+NNN)/2)).*fil(1:floor((2+NNN)/2));
for i=floor((2+NNN)/2):NNN
    invf(i)=conj(invf(NNN+2-i));
end
po=po-real(ifft(invf));
%bandpass filter, only process the noisy frequency part po

pp1=reshape(po,n,floor(NNN/n));%break data into serveral periods
pp2=0*po';%store the estimated sinusoidal noise

for i=0:lr/2
    pt1=pp1(:,1:lr);
    if rd==1
        c=srt(n,lr);
        p11=pt1*0;
        for k=1:n
            p11(k,:)=pt1(k,c(k,:));
        end
        pt1=p11;
    end%randomization
    [u,s,v]=svd(pt1);
    ss=s*0;ss(1:rk,1:rk)=s(1:rk,1:rk)*rto;%svd
    pt2=u*ss*v';
    xx=sum(pt2,2)/lr;%average
    pp2(i*n+1:(i+1)*n)=xx;
end
%begining periods


for i=lr/2+1:floor(NN/n)-lr/2
    pt1=pp1(:,i-lr/2:i+lr/2);
    if rd==1
        c=srt(n,lr+1);
        p11=pt1*0;
        for k=1:n
            p11(k,:)=pt1(k,c(k,:));
        end
        pt1=p11;
    end
    [u,s,v]=svd(pt1);
    ss=s*0;ss(1:rk,1:rk)=s(1:rk,1:rk)*rto;
    pt2=u*ss*v';
    xx=sum(pt2,2)/(lr+1);%
    pp2(i*n+1:(i+1)*n)=xx;
end
%middle periods

for i=floor(NN/n)-lr/2+1:floor(NN/n)
    pt1=pp1(:,floor(NN/n)-lr+1:floor(NN/n));
    if rd==1
        c=srt(n,lr);
        p11=pt1*0;
        for k=1:n
            p11(k,:)=pt1(k,c(k,:));
        end
        pt1=p11;
    end
    [u,s,v]=svd(pt1);
    ss=s*0;ss(1:rk,1:rk)=s(1:rk,1:rk)*rto;
    pt2=u*ss*v';
    xx=sum(pt2,2)/lr;%
    pp2(i*n+1:(i+1)*n)=xx;
end
%end periods

x2=po-pp2';%subtract the sinusoidal noise 
p4=x2+real(ifft(invf));
%add the other part beyond the bandpass range back

p5=p4(lr/2*n+1:end-lr/2*n);
p5=p5(1:NN-lr*n);
p4=p5(1:ninter:end);
%remove the effects of interpolation and padding
end

