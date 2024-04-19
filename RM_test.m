clc;
clear all; %#ok<CLALL>
close all;
%%
frameCount = 1;
ID=4; %ncell_ID
mu = 0;
ss=SsGenerator();
[pss,sss]=ss.getSsSignalsByCellInfo(ID);
r=ResourceMapper();
r.createResourceGrid(mu,frameCount);

r.addSsBlockByCase('A',[0,1],ID,pss,sss,ones(100,432)*20,ones(100,144),0,0);

h=heatmap(abs(r.resourceGrid(:,1:50)));
grid off
