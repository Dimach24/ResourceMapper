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

r.addSsBlockToResourceGrid(ID,pss,sss,ones(1,432)*20,ones(1,144),5,0);

h=heatmap(abs(r.resourceGrid(:,1:50)));
grid off
