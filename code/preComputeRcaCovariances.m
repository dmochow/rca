function [sumXX,sumYY,sumXY,nPointsInXX,nPointsInYY,nPointsInXY]=preComputeRcaCovariances(data,condRange,subjRange)
% [SUMXX,SUMYY,SUMXY,NPOINTSINXX,NPOINTSINYY,NPOINTSINXY]=PRECOMPUTERCACOVARIANCES(DATA)
% 
% take in epoched data volume(s) to produce within- and
% across-trial covariance matrices for input into rcaTrain()
%
% data can either be a 3-D array (samples x channels x trials), or a cell
% array (conditions x subjects) of 3-D data arrays (as above)
%
% throughout data, NaNs indicate missing data values (i.e., artifacts)
%
% sumXX, sumYY, sumXY: unnormalized auto- and cross-covariances
% nPointsInXX, nPointsInYY, nPointsInXY: number of non-missing data points used to compute covariances
%
% covariances are left unnormalized to allow for subsequent averaging
% across subjects or conditions, as appropriate
%
% (c) Jacek P. Dmochowski, 2014
if ~iscell(data), 
    nCond=1; nSubjects=1; subjRange=1; condRange=1; data={data}; 
else
    if nargin<3, subjRange=1:size(data,2); end
    if nargin<2, condRange=1:size(data,1); end
    data=data(condRange,subjRange);
    nCond=size(data,1);
    nSubjects=size(data,2);
end

fprintf('Selected %d subjects and %d conditions for training... \n',nSubjects,nCond);

[~,nElectrodes,~]=size(data{1,1});  % assume uniform (?)
sumXX=zeros(nCond,nSubjects,nElectrodes,nElectrodes);sumYY=zeros(nCond,nSubjects,nElectrodes,nElectrodes);sumXY=zeros(nCond,nSubjects,nElectrodes,nElectrodes);
nPointsInXX=zeros(nCond,nSubjects,nElectrodes,nElectrodes);nPointsInYY=zeros(nCond,nSubjects,nElectrodes,nElectrodes);nPointsInXY=zeros(nCond,nSubjects,nElectrodes,nElectrodes);

for cond=1:nCond
    for subj=1:nSubjects
        fprintf('Computing covariances for subject %d and condition %d... \n',subjRange(subj),condRange(cond));
        thisVolume=data{cond,subj};
        [nSamples,nElectrodes,~]=size(thisVolume);
        if nSamples<nElectrodes, warning('Number of samples is less than the number of electrodes'); end
        
        
        nTrials=size(thisVolume,3);
        pindx=combnk(1:nTrials,2);
        pindx= cat(1,pindx,pindx(:,[2 1])); % Lucas' trick to ensure that Rxx=Ryy
        
        % offset by mean local to this subject-condition
        concatX=permute(thisVolume(:,:,pindx(:,1)),[2 1 3]);
        concatY=permute(thisVolume(:,:,pindx(:,2)),[2 1 3]);
        concatX=concatX(:,:); concatY=concatY(:,:);
        concatX=concatX-repmat(nanmean(concatX,2),[1 size(concatX,2)]);
        concatY=concatY-repmat(nanmean(concatY,2),[1 size(concatY,2)]);
        
        nPointsInXX(cond,subj,:,:)= double(~isnan(concatX))*double(~isnan(concatX))';
        nPointsInYY(cond,subj,:,:)= double(~isnan(concatY))*double(~isnan(concatY))';
        nPointsInXY(cond,subj,:,:)= double(~isnan(concatX))*double(~isnan(concatY))';
        
        concatX(isnan(concatX))=0; concatY(isnan(concatY))=0;
        
        sumXX(cond,subj,:,:)=concatX*concatX'; sumYY(cond,subj,:,:)=concatY*concatY'; sumXY(cond,subj,:,:)=concatX*concatY';
    end
end

sumXX=squeeze(sumXX); sumYY=squeeze(sumYY); sumXY=squeeze(sumXY);
nPointsInXX=squeeze(nPointsInXX); nPointsInYY=squeeze(nPointsInYY);  nPointsInXY=squeeze(nPointsInXY);

