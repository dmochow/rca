clear all; close all; clc
addpath(genpath('.'))
load ./data/rcaDemoData

%% run with default parameters
[dataOut,W,A] = rcaRun(myCellData);

%% run with custom parameters
%[dataOut,W,A] = rcaRun(myCellData,7,5,[1 3],[1 3 4],1);  % diagonalize 7 dimensions, compute 5 components, train on subset of conditions/subjects

%% run for one subject, all conditions
%[dataOut,W,A] = rcaRun(myCellData(:,1));

%% run for all subjects, one condition
%[dataOut,W,A] = rcaRun(myCellData(1,:));

%% run for one subject, one condition
%[dataOut,W,A] = rcaRun(myCellData{1,1});
%[dataOut,W,A] = rcaRun(myCellData{1,1},7,5,[],[],1);
%[dataOut,W,A] = rcaRun(myCellData{1,1},10,[],[],[],1);