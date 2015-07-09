clear all
addpath(genpath('.'))
load('./data/rcaDemoData','data','condLeft','condRight','fs');  

data

%% run with default parameters
[dataOut,W,A] = rcaRun(data);

%% make a plot of RCA results
grpIndx={[1 4],[2 3]};
condLabels={'easy','hard'};
colors={'b','r'};
tvec=1:2*fs;  % 2-second trials
[h,hlg,grpData,grpDataSem] = rcaPlot( dataOut , A , grpIndx  , tvec , colors, condLabels ); % plot results

% %% run with custom parameters
% nReg=5;
% nComp=4;
% [dataOut,W,A] = rcaRun(data,nReg,nComp);  % diagonalize 5 dimensions, compute 4 components
% 
% %% learn only on the "left conditions"
% nReg=5; 
% nComp=3;
% [dataOutLeft,Wleft,Aleft] = rcaRun(data,nReg,nComp,condLeft);
% 
% %% learn on "right conditions" and then apply resulting weights to "left conditions"
% nReg=5; 
% nComp=3;
% [dataOutRight,Wright,Aright] = rcaRun(data,nReg,nComp,condRight);
% dataOutLeftRight = rcaProject(data(condLeft),Wright);  % apply to conditions 1 and 3
% [h,hlg,grpData,grpDataSem] = rcaPlot( dataOutLeftRight , Aright , {1,2}  , 1:2*fs , {'b','r'}, {'easy','hard'} );
