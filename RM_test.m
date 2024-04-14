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


r.addPbchDmRsToResourceGrid(ID,ones(frameCount*2^mu*10,144)*40);
r.addPbchToResourceGrid(ID,ones(frameCount*2^mu*10,432)*20);
r.addPssToResourceGrid(pss);
r.addSssToResourceGrid(sss);
h=heatmap(abs(r.resourceGrid(:,1:50)));
grid off
labels=string(1:240);
labels(mod(1:240,25)~=0)='';
h.YDisplayLabels=labels;
h.YDisplayData=flipud(h.YDisplayData);