%% demo for the AFNLM filtering
clear;close all;
%% load noisy and clean data
load('a.mat');
load('a1.mat');% noisy
%% AFNLM filter

fmi=1;fmx=70;%frequency bandwidth in sample-points unit
l1=100;l2=150;%filtering parameter
f=10;s=5;%size of searching area and patch
sigma=9000;

tic;decon = fx_decon(a1,0.001,40,0.01,1,200);toc;%pre-processed data
deconnoi=a1-decon;


[x,noi]=afnlm(a1,f,s,l1,l2,fmi,fmx,sigma,decon,deconnoi);

%% plot and SNR
lim1=-1;lim2=1;
x1=100;y1=100;dx=700;dy=550;

figure;imagesc(a);ax = gca;colormap(gray);set(ax, 'CLim', [lim1 lim2]);set(gcf,'position',[x1,y1,dx,dy]);colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);
figure;imagesc(a1);ax = gca;colormap(gray);set(ax, 'CLim', [lim1 lim2]);set(gcf,'position',[x1,y1,dx,dy]);colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);
figure;imagesc(x);ax = gca;colormap(gray);set(ax, 'CLim', [lim1 lim2]);set(gcf,'position',[x1,y1,dx,dy]);colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);
figure;imagesc(decon);ax = gca;colormap(gray);set(ax, 'CLim', [lim1 lim2]);set(gcf,'position',[x1,y1,dx,dy]);colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);
figure;imagesc(noi);ax = gca;colormap(gray);set(ax, 'CLim', [lim1 lim2]);set(gcf,'position',[x1,y1,dx,dy]);colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);
figure;imagesc(deconnoi);ax = gca;colormap(gray);set(ax, 'CLim', [lim1 lim2]);set(gcf,'position',[x1,y1,dx,dy]);colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

snro=SNR(a,a1);
snrxaf=SNR(a,x);
snrd=SNR(a,decon);