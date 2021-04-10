clear; close all;
%% load synthetic data
load F_s1.mat;
load F_sn1.mat;
s=s';% clean data
s3=s3';% noise
sn=s+s3;% noisy data

x1=200;y1=200;dx=1500;dy=400;
figure;
plot(1:length(s),[s,sn],'LineWidth',1.5);ax=gca;set(gcf,'position',[x1,y1,dx,dy]);
xlabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);xlim([0,length(s)])
ylabel('Amplitude','FontName','Arial','FontWeight','Bold','FontSize',14);ylim([-6,6])
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1.5);
% xticks(0:100:length(s));
% xticklabels({'0','238.1','476.19','714.29','952.38','1190.5'});
legend({'noise-free trace','noisy trace'},'FontName',...
    'Arial','FontSize',16,'LineWidth',2,'Location','southwest');

c=abs(fft([s,sn]));
d=angle(fft([s,sn]));
figure;plot(1:100,log(c(1:100,:)),'LineWidth',1.5);ax=gca;set(gcf,'position',[x1,y1,dx,dy]);
xlabel('Frequency (Hz)','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Log amplitude','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1.5);
legend({'noise-free amplitude spectrum','noisy amplitude spectrum'},'FontName',...
    'Arial','FontSize',16,'LineWidth',2,'Location','southeast');

figure;plot(1:100,d(1:100,:),'LineWidth',1.5);ax=gca;set(gcf,'position',[x1,y1,dx,dy]);
xlabel('Frequency (Hz)','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Phase','FontName','Arial','FontWeight','Bold','FontSize',14);ylim([-4,4]);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1.5);
legend({'noise-free phase spectrum','noisy phase spectrum'},'FontName',...
    'Arial','FontSize',16,'LineWidth',2,'Location','southwest');
%% denoise

ninter=21;%the number of the interpolation points
dt=0.001;% time interval
lr=26;% the number of periods used for the RPCA matrix (lateral size of C)
rk=3;% truncted rank value of C
wid=40;%bandpass filter radius in sample-points unit
rto=0.98;% scaling factor for the ranks
rd=0;% if using randomization 0=no, 1=yes
[p4,ffs]= func_rpca(sn,ninter,lr,rk,wid,rto,rd,dt);


x1=200;y1=200;dx=1500;dy=400;
figure;
hold on;plot(p4,'LineWidth',1.5);plot(s,'r','LineWidth',1.5);set(gcf,'position',[x1,y1,dx,dy]);
xlabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Amplitude','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1.5,'Ylim',[-1.5,1.5]);
legend({'denoised result of RPCA','noise-free trace'},'FontName',...
    'Arial','FontSize',16,'LineWidth',2,'Location','southwest');

figure;
plot(p4-s,'LineWidth',1.5);set(gcf,'position',[x1,y1,dx,dy]);
xlabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Amplitude','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1.5,'Ylim',[-1.5,1.5]);
legend({'error between estimation of RPCA and real signal'},'FontName',...
    'Arial','FontSize',16,'LineWidth',2,'Location','southwest');

ccc=abs(fft([p4,s]));
ddd=angle(fft([p4,s]));
figure;plot(1:100,log(ccc(1:100,:)),'LineWidth',1.5);set(gcf,'position',[x1,y1,dx,dy]);
xlabel('Frequency (Hz)','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Log amplitude','FontName','Arial','FontWeight','Bold','FontSize',14);ylim([-2,5]);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1.5);
legend({'amplitude spectrum of RPCA','noise-free amplitude spectrum'},'FontName',...
    'Arial','FontSize',16,'LineWidth',2,'Location','southeast');

figure;plot(1:100,ddd(1:100,:),'LineWidth',1.5);set(gcf,'position',[x1,y1,dx,dy]);
xlabel('Frequency (Hz)','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Phase','FontName','Arial','FontWeight','Bold','FontSize',14);ylim([-4,4]);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1.5);
legend({'phase spectrum of RPCA','noise-free phase spectrum'},'FontName',...
    'Arial','FontSize',16,'LineWidth',2,'Location','southwest');

snro=snr(s,sn-s);
snrr=snr(s,p4-s);