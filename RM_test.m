clc;
clear all; %#ok<CLALL>
close all;
addpath(genpath(pwd));
%%
frameCount = 6;
ID=4; %ncell_ID
mu = 1;
[pss,sss]=SsGenerator.getSsSignalsByCellInfo(ID);
r=ResourceMapper();
r.createResourceGrid(mu,frameCount,false,30);

r.addSsBlockByCase('B',0:8,ID,pss,sss,ones(5*500,432)*20,ones(5*500,144)*40,0,0);

s=pcolor(abs(r.resourceGrid(1:260,1:80)));
s.EdgeColor='none';
grid off
a=gca();
a.YDir='normal';
xlabel('l (номер OFDM символа)')
ylabel('k (номер поднесущей)')