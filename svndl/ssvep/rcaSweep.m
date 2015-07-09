function [rcaDataReal,rcaDataImag,binData]=rcaSweep(pathnames,binsToUse,freqsToUse,nReg,nComp,show)
% perform RCA on sweep SSVEP data exported to RLS or DFT format
%
% [RCADATAREAL,RCADATAIMAG,BINDATA]=RCASWEEP(PATHNAMES,[BINSTOUSE],[FREQSTOUSE],[NREG],[NCOMP],[SHOW])
%
% pathnames (required): cell vector of string directory names housins DFT_c00x.txt or RLS_c00x.txt exports
%   for example,  pathnames={'/Volumes/Denali_4D2/rca/s001/','/Volumes/Denali_4D2/rca/s002/','/Volumes/Denali_4D2/rca/s003/'}
% binsToUse: vector of bin indices to include in RCA (defaults to bin 0 or average across bins)
% freqsToUse: vector of frequency indices to include in RCA (defaults to 1 or first harmonic ?)
% nReg: RCA regularization parameter (defaults to 9)
% nComp: number of RCs to retain (defaults to 3)
% show: 1 to see a figure of sweep amplitudes for each harmonic and component (defaults to 1), 0 to not display
%
% rcaDataReal and rcaDataImag are cell arrays of size
% nConditions-by-nSubjects, where each element of the cell array has
% dimensions coefficient-by-component-by-trial
% NB: the coefficients are arranged first across bins, and then
% concatenated across frequencies.  For example, if one has 10 bins and 2 frequencies, than the
% 20-element coefficient vector refers to
% [(bin1,f1),(bin2,f1),...,(bin10,f1),(bin1,f2),...(bin10,f2)]
%
% bindata is a struct containing subject and condition averaged data that
% contains the following fields:
%
%   realBins: bin-by-harmonic-by-component array of REAL RC coefficients
%   imagBins: bin-by-harmonic-by-component array of IMAGINARY RC coefficients
%   ampBins: bin-by-harmonic-by-component array of RC amplitudes
%   phaseBins: bin-by-harmonic-by-component array of RC phases
%
%   noise1RealBins: bin-by-harmonic-by-component array of noise-band (lo) real RC coefficients
%   noise1ImagBins: bin-by-harmonic-by-component array of noise-band (lo) imaginary RC coefficients
%   noise1AmpBins: bin-by-harmonic-by-component array of noise-band (lo) RC amplitudes
%   noise1PhaseBins: bin-by-harmonic-by-component array of noise-band (lo) RC phases
%
%   noise2RealBins: bin-by-harmonic-by-component array of noise-band (hi) real RC coefficients
%   noise2ImagBins: bin-by-harmonic-by-component array of noise-band (hi) imaginary RC coefficients
%   noise2AmpBins: bin-by-harmonic-by-component array of noise-band (hi) RC amplitudes
%   noise2PhaseBins: bin-by-harmonic-by-component array of noise-band (hi) RC phases
%
% Jacek P. Dmochowski, 2015, report bugs to dmochowski@gmail.com


if nargin<6, show=1; end
if nargin<5, nComp=3; end
if nargin<4, nReg=9; end
if nargin<3, freqsToUse=1; end
if nargin<2, binsToUse=0; end
if nargin<1, error('Must specify at least one input argument'); end

%% if pathanmes is a string, convert to cell
if ~iscell(pathnames)
    try
        pathnames=mat2cell(pathnames);
    catch
        error('unable to parse pathnames: check that it is a cell array of strings');
    end
end
%% read in signal and noise data
nSubjects=numel(pathnames);
sensorData={};
noiseData1={};
noiseData2={};
freqIndices=cell(nSubjects,1);
binIndices=cell(nSubjects,1);
for s=1:nSubjects
    [signalData,indF,indB,noise1,noise2]=textExportToRca(pathnames{s},binsToUse,freqsToUse);
    freqIndices{s}=indF;
    binIndices{s}=indB;
    sensorData(:,s)=signalData;
    cellNoiseData1(:,s)=noise1;
    cellNoiseData2(:,s)=noise2;
end

%% check freq and bin indices for consistency across subjects
for s=1:nSubjects
    if sum(abs(freqIndices{s}-freqIndices{1}))~=0 && sum(abs(binIndices{s}-binIndices{1}))~=0
        error('Frequency and bins indices vary across subjects: check consistency of DFT/RLS exports');
    end
end

% if they're all the same, make into a regular array
freqIndices=freqIndices{1};
binIndices=binIndices{1};

%% run RCA
[rcaData,W,A]=rcaRun(sensorData,nReg,nComp);
noiseData1=rcaProject(cellNoiseData1,W);
noiseData2=rcaProject(cellNoiseData2,W);

% create a "component" of just channel Oz for performance evaluation
wOz=zeros(128,1); wOz(75)=1; % this is only valid with the EGI HydroCel 128 System
ozData=rcaProject(sensorData,wOz);
ozNoiseData1=rcaProject(cellNoiseData1,wOz);
ozNoiseData2=rcaProject(cellNoiseData2,wOz);

%% convert to real/imaginary
rcaDataReal=cell(size(rcaData));
rcaDataImag=cell(size(rcaData));

noiseData1Real=cell(size(noiseData1));
noiseData1Imag=cell(size(noiseData1));

noiseData2Real=cell(size(noiseData2));
noiseData2Imag=cell(size(noiseData2));

ozDataReal=cell(size(ozData));
ozDataImag=cell(size(ozData));

ozNoiseData1Real=cell(size(ozNoiseData1));
ozNoiseData1Imag=cell(size(ozNoiseData1));

ozNoiseData2Real=cell(size(ozNoiseData2));
ozNoiseData2Imag=cell(size(ozNoiseData2));

% loop through elements of cell array and cut into real and imag halves
for i=1:size(rcaDataReal,1)
    for j=1:size(rcaDataReal,2)
   
        % check number of coefficients in this cell
        thisData=rcaData{i,j};
        nSamples=size(thisData,1);
        if rem(nSamples,2)~=0, error('odd number of samples is not consistent with ssvep data');  end
        
        % cut and store
        rcaDataReal{i,j}=rcaData{i,j}(1:nSamples/2,:,:);
        rcaDataImag{i,j}=rcaData{i,j}(nSamples/2+1:nSamples,:,:);
        
        % assumes noiseData1 and noiseData2 have same size as rcaData
        noiseData1Real{i,j}=noiseData1{i,j}(1:nSamples/2,:,:);
        noiseData1Imag{i,j}=noiseData1{i,j}(nSamples/2+1:nSamples,:,:);
        
        noiseData2Real{i,j}=noiseData2{i,j}(1:nSamples/2,:,:);
        noiseData2Imag{i,j}=noiseData2{i,j}(nSamples/2+1:nSamples,:,:);
        
        ozDataReal{i,j}=ozData{i,j}(1:nSamples/2,:,:);
        ozDataImag{i,j}=ozData{i,j}(nSamples/2+1:nSamples,:,:);
        
        ozNoiseData1Real{i,j}=ozNoiseData1{i,j}(1:nSamples/2,:,:);
        ozNoiseData1Imag{i,j}=ozNoiseData1{i,j}(nSamples/2+1:nSamples,:,:);
        
        ozNoiseData2Real{i,j}=ozNoiseData2{i,j}(1:nSamples/2,:,:);
        ozNoiseData2Imag{i,j}=ozNoiseData2{i,j}(nSamples/2+1:nSamples,:,:);
        
    end
end


% from this point on, the code tries to compile some basic features and
% display a few curves

% in principle, all the data that one needs for their analysis is in
% rcaDataReal and rcaDataImag, and perhaps
% noiseData1Real/noiseData1Imag/noiseData2Real/noiseData2Imag

%% aggregate across subjects AND conditions
rcaDataRealAll=cat(3,rcaDataReal{:});
rcaDataImagAll=cat(3,rcaDataImag{:});

noiseData1RealAll=cat(3,noiseData1Real{:});
noiseData1ImagAll=cat(3,noiseData1Imag{:});

noiseData2RealAll=cat(3,noiseData2Real{:});
noiseData2ImagAll=cat(3,noiseData2Imag{:});

ozDataRealAll=cat(3,ozDataReal{:});
ozDataImagAll=cat(3,ozDataImag{:});

ozNoiseData1RealAll=cat(3,ozNoiseData1Real{:});
ozNoiseData1ImagAll=cat(3,ozNoiseData1Imag{:});

ozNoiseData2RealAll=cat(3,ozNoiseData2Real{:});
ozNoiseData2ImagAll=cat(3,ozNoiseData2Imag{:});

%% trial average for sweep plots
muRcaDataReal=nanmean(rcaDataRealAll,3);
muRcaDataImag=nanmean(rcaDataImagAll,3);

muNoiseData1Real=nanmean(noiseData1RealAll,3);
muNoiseData1Imag=nanmean(noiseData1ImagAll,3);

muNoiseData2Real=nanmean(noiseData2RealAll,3);
muNoiseData2Imag=nanmean(noiseData2ImagAll,3);

muOzDataReal=nanmean(ozDataRealAll,3);
muOzDataImag=nanmean(ozDataImagAll,3);

muOzNoiseData1Real=nanmean(ozNoiseData1RealAll,3);
muOzNoiseData1Imag=nanmean(ozNoiseData1ImagAll,3);

muOzNoiseData2Real=nanmean(ozNoiseData2RealAll,3);
muOzNoiseData2Imag=nanmean(ozNoiseData2ImagAll,3);

%% make a frequency by bin by trial cube
nFreqs=numel(unique(freqIndices));
nBins=numel(unique(binIndices));
realBins=zeros(nBins,nFreqs,nComp);
imagBins=zeros(nBins,nFreqs,nComp);
noise1RealBins=zeros(nBins,nFreqs,nComp);
noise1ImagBins=zeros(nBins,nFreqs,nComp);
noise2RealBins=zeros(nBins,nFreqs,nComp);
noise2ImagBins=zeros(nBins,nFreqs,nComp);
ozRealBins=zeros(nBins,nFreqs,nComp);
ozImagBins=zeros(nBins,nFreqs,nComp);
ozNoise1RealBins=zeros(nBins,nFreqs,nComp);
ozNoise1ImagBins=zeros(nBins,nFreqs,nComp);
ozNoise2RealBins=zeros(nBins,nFreqs,nComp);
ozNoise2ImagBins=zeros(nBins,nFreqs,nComp);
for c=1:nComp
    for f=1:nFreqs
        realBins(:,f,c)=muRcaDataReal(freqIndices==freqsToUse(f),c);
        imagBins(:,f,c)=muRcaDataImag(freqIndices==freqsToUse(f),c);
        
        noise1RealBins(:,f,c)=muNoiseData1Real(freqIndices==freqsToUse(f),c);
        noise1ImagBins(:,f,c)=muNoiseData1Imag(freqIndices==freqsToUse(f),c);
        
        noise2RealBins(:,f,c)=muNoiseData2Real(freqIndices==freqsToUse(f),c);
        noise2ImagBins(:,f,c)=muNoiseData2Imag(freqIndices==freqsToUse(f),c);
        
        ozRealBins(:,f,c)=muOzDataReal(freqIndices==freqsToUse(f),1); % just repeat for Oz
        ozImagBins(:,f,c)=muOzDataImag(freqIndices==freqsToUse(f),1);
        
        ozNoise1RealBins(:,f,c)=muOzNoiseData1Real(freqIndices==freqsToUse(f),1);
        ozNoise1ImagBins(:,f,c)=muOzNoiseData1Imag(freqIndices==freqsToUse(f),1);
        
        ozNoise2RealBins(:,f,c)=muOzNoiseData2Real(freqIndices==freqsToUse(f),1);
        ozNoise2ImagBins(:,f,c)=muOzNoiseData2Imag(freqIndices==freqsToUse(f),1);
        
    end
end

ampBins=sqrt(realBins.^2+imagBins.^2);
phaseBins=atan(imagBins./realBins);

noise1AmpBins=sqrt(noise1RealBins.^2+noise1ImagBins.^2);
noise1PhaseBins=atan(noise1ImagBins./noise1RealBins);

noise2AmpBins=sqrt(noise2RealBins.^2+noise2ImagBins.^2);
noise2PhaseBins=atan(noise2ImagBins./noise2RealBins);

ozAmpBins=sqrt(ozRealBins.^2+ozImagBins.^2);
ozPhaseBins=atan(ozImagBins./ozRealBins);

ozNoise1AmpBins=sqrt(ozNoise1RealBins.^2+ozNoise1ImagBins.^2);
ozNoise1PhaseBins=atan(ozNoise1ImagBins./ozNoise1RealBins);

ozNoise2AmpBins=sqrt(ozNoise2RealBins.^2+ozNoise2ImagBins.^2);
ozNoise2PhaseBins=atan(ozNoise2ImagBins./ozNoise2RealBins);

snrBins=2*ampBins./(noise1AmpBins+noise2AmpBins);
ozSnrBins=2*ozAmpBins./(ozNoise1AmpBins+ozNoise2AmpBins);

binData.realBins=realBins;
binData.imagBins=imagBins;
binData.ampBins=ampBins;
binData.phaseBins=phaseBins;

binData.noise1RealBins=noise1RealBins;
binData.noise1ImagBins=noise1ImagBins;
binData.noise1AmpBins=noise1AmpBins;
binData.noise1PhaseBins=noise1PhaseBins;

binData.noise2RealBins=noise2RealBins;
binData.noise2ImagBins=noise2ImagBins;
binData.noise2AmpBins=noise2AmpBins;
binData.noise2PhaseBins=noise2PhaseBins;

binData.ozRealBins=ozRealBins;
binData.ozImagBins=ozImagBins;
binData.ozAmpBins=ozAmpBins;
binData.ozPhaseBins=ozPhaseBins;

binData.ozNoise1RealBins=ozNoise1RealBins;
binData.ozNoise1ImagBins=ozNoise1ImagBins;
binData.ozNoise1AmpBins=ozNoise1AmpBins;
binData.ozNoise1PhaseBins=ozNoise1PhaseBins;

binData.ozNoise2RealBins=ozNoise2RealBins;
binData.ozNoise2ImagBins=ozNoise2ImagBins;
binData.ozNoise2AmpBins=ozNoise2AmpBins;
binData.ozNoise2PhaseBins=ozNoise2PhaseBins;

binData.snrBins=snrBins;
binData.ozSnrBins=ozSnrBins;


%%

if show
    figure
    for f=1:nFreqs
        for c=1:nComp
            %subplot(nFreqs,nComp,(f-1)*nComp+c); hold on
            subplot(nFreqs,nComp+1,(f-1)*(nComp+1)+c); hold on
            plot(ampBins(:,f,c),'--ok','MarkerFaceColor','k');
            plot(noise1AmpBins(:,f,c),'--or','MarkerFaceColor','r');
            plot(noise2AmpBins(:,f,c),'--og','MarkerFaceColor','g');
            
            if f==1, title(['Component ' num2str(c)']); end
            if c==1, ylabel(['Frequency ' num2str(freqsToUse(f))']); end
            if f==1 && c==1, hlg=legend('signal','noise1','noise2'); end
            
            subplot(nFreqs,nComp+1,f*(nComp+1)); hold on
            plot(ozAmpBins(:,f,c),'--ok','MarkerFaceColor','k');
            plot(ozNoise1AmpBins(:,f,c),'--or','MarkerFaceColor','r');
            plot(ozNoise2AmpBins(:,f,c),'--og','MarkerFaceColor','g');
            if f==1, title('Electrode Oz'); end
          
            
        end
    end
    set(hlg,'box','off');
    
    figure % SNR comparison
    for f=1:nFreqs     
        subplot(nFreqs,1,f); hold on
        plot(snrBins(:,f,1),'--ok','MarkerFaceColor','k');
        plot(ozSnrBins(:,f,1),'--or','MarkerFaceColor','r');
        title(['Frequency ' num2str(freqsToUse(f))']);
        ylabel('SNR');
        if f==1, hlg=legend('RC1','Oz'); end
    end
    set(hlg,'box','off');
    
    
    
end