function ffs=precise_dft(x,f1,f2,n,ninter,dt)
df=(f2-f1)/n;
dt=dt/ninter;
m=length(x);
ff=zeros(n,1);
for j=0:n-1
et=exp(-1i*2*pi*(df*j+f1).*(1:m)*dt);
ff(j+1)=sum(x'.*et);
end

ffabs=abs(ff);
ffs=(find(ffabs==max(ffabs))-1)*df+f1;%
end