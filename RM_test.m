clc;
clear all; %#ok<CLALL>
close all;
%%
frameCount = 6;
ID=4; %ncell_ID
mu = 4;
[pss,sss]=SsGenerator.getSsSignalsByCellInfo(ID);
r=ResourceMapper();
r.createResourceGrid(mu,frameCount);

r.addSsBlockByCase('G',0:8,ID,pss,sss,ones(5*500,432)*20,ones(5*500,144)*40,0,0);

h=heatmap(abs(r.resourceGrid(1:300,100:end)));
grid off
