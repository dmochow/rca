clear all; close all; clc;
addpath(genpath('/Users/jacek/Dropbox/rca'));

%% take powerdiva export and convert it to cell array format
dataPath='/Users/jacek/Documents/MATLAB/EEG/dirDisc/data/MT/';
%dataPath='/Users/jacek/Documents/MATLAB/EEG/appMotion/data/S1/';
%dataPath='/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/16_groups/Session_1/20140212_PJK/Kohler_20140212_1522/Exp_MATL_HCN_128_Avg/';
cellData=exportToRcaReady(dataPath,1,1,[20 40]);


%% %% run RCA
nReg=7;
nComp=5;
[rcaData,W,A]=rcaRun(cellData,nReg,nComp);










