function [new,tt]=nlm_f(a,x,noi,f,s,vara,l1,l2,fmi,fmx)
% nlm in frequency slice
% INPUT
% a: input noisy data
% f:  the half side length of the searching area
% s:  the half side length of the patch area
% l1,l2: filtering parameter
% fmi,fmx: frequency bandwidth in sample points
% vara: variance of gauss kernal
% x: denoised data by the pre-proecssing
% noi:  noise data by the pre-proecssing

% OUTPUT
% new:  denoised data
% tt: adaptive weights
%
% Copyright (C) 2018 Hang Wang
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%  Reference:   
%
%  [1]Buades, A., B. Coll, and J. M. Morel, 2005, A review of image denoising
%  algorithms, with a new one: Multiscale Modeling and Simulation


%the size of input data and denoising data
[m,n]=size(a);
new=zeros(m,n);

%padding the original data and get the range of pixels that wait to be calculated
a1=padarray(a,[0,f],'symmetric');
m1=1;
m2=m;
n1=f+1;
n2=n+f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive h
xfabs=mean(abs(x(fmi:fmx,:)),2);
noifabs=mean(abs(noi(fmi:fmx,:)),2);
xfabs=repmat(xfabs,1,n+f);
noifabs=repmat(noifabs,1,n+f);
tt=(noifabs+0.001)./(xfabs+0.001);
zmax=max(max(tt));
zmin=min(min(tt));
h1=(l2-l1)*(tt-zmin)/(zmax-zmin)+l1;
for i=1:size(h1,2)
    h1(:,i)=smooth(h1(:,i),10);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=m1:m2%loop of slices
    for j=n1:n2%loop of searching windows
        
        %searching window
        ftemp=a1(i,j-f:j+f);
        
        %reference patch
        refer=ftemp(f+1-s:f+1+s);
        [~,q]=size(ftemp);
        y1=s+1;
        y2=q-s;
        d=zeros(1,q-2*s);
        
        % loop of patches
            for k2=y1:y2
                stemp=ftemp(k2-s:k2+s);
                % calculate the distance between reference patch and the choosed patch and gauss kernal
                [d(k2-s),gs]=distan(refer,stemp,vara);
            end
        
        %normalize the weight and get the estimation of this pixel
        Z=sum(sum(exp(-d./(h1(i,j)^2))));
        w=1/Z.*exp(-d./(h1(i,j)^2));
        sigm=sum(w.*ftemp(y1:y2));
        new(i,j-f)=sigm;
    end
    disp(i)
end
end

%calculate the distance of patches in different location
function [dw,gs]=distan(a,b,devi)

%initialization of input patches
[~,n]=size(a);
gs=a*0;
c=a-b;

%generate the gauss kernal
cent2=floor(n/2);

    for j=1:n
        gs(j)=exp(-((j-cent2).^2)./2/devi);
    end

%calculate the distance under the gauss kernal
dw=(norm(gs.*c,2))^2;
end