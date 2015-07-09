function [cellData,timelines]=exportToRcaReady(dataPath,removeEyes,nanArtifacts,notchFreqs)
% CELLDATA=EXPORTTORCAREADY(DATAPATH,[REMOVEEYES],[NANARTIFACTS],[NOTCHFREQS])
%
% This function converts a PowerDiva Matlab export into a data structure
% that can be analyzed with RCA.
%
% dataPath: pathname to the PowerDiva export (must have RtSeg and
% raw_cxx_tyy.mat files)
%
% removeEyes:  1 (default) to regress out EOG channels from EEG, 0 to not
% regress out
%
% nanArtifacts: 1 (default) to insert NaNs into artifactual samples
% identified by PowerDiva; 0 to not insert NaNs
%
% notchFreqs: a vector of frequencies (in Hz) that are to be notched out of the data (defaults to [])
% 
% Jacek P. Dmochowski, Peter J. Kohler, 2015
%
% contact dmochowski@gmail.com for support

%% sweep-style, fixed-style, erp-style

%% spoof-Axx 


%% separate function to trim time window of cell array
if nargin<4, notchFreqs=[]; end;
if nargin<3, nanArtifacts=1; end;
if nargin<2, removeEyes=1; end;

bslIndx=1:20; % hard-coded (DWIB)

curDir=pwd; % save for later (return to original path after routine)

cd(dataPath);

RTsegFiles = subfiles('RTSeg_*.mat');
nSeg=numel(RTsegFiles);
condInds=cell(nSeg,1);
segCount=0;

for z=1:nSeg
    load(RTsegFiles{z}); % load data
    if size(TimeLine,1)>0  % ignore aborted runs
        segCount = segCount+1;
        
        try
            timelines{segCount}=TimeLine;
        catch
            warning('unable to find timeline for this segment');
        end
        
        condInds{segCount}=unique([TimeLine.cndNmb]);  % store condition indices for this run
        
        
        rawdata=cell(numel(condInds{segCount}),1);
        for t=1:size(TimeLine,1)
            rawFileName=['Raw_c' num2str(TimeLine(t).cndNmb,'%03.0f') '_t' num2str(TimeLine(t).trlNmb,'%03.0f')];
            load(rawFileName);
            rawtrial= double(RawTrial).*repmat(Ampl.',size(RawTrial,1),1) +   repmat( Shift.' , size(RawTrial,1) , 1 );  % convert to volts
            
            
            if removeEyes
                X=rawtrial(:,1:end-2); % data channels
                V=rawtrial(:,end-1:end);  % HEOG/VEOG
                if ~isempty(find(isnan(V), 1))
                    warning('EOG channels contain NaNs')
                else
                    rawtrial=( eye(size(X,1))-V*pinv(V) )*X; %A=pinv(V)*X; % transfer function from eyes to data electrodes
                end
            else
                rawtrial=rawtrial(:,1:end-2);  % discard EOG channels
            end
            
            if ~isempty(notchFreqs)
                for f=1:numel(notchFreqs)

                    w0=notchFreqs(f)/(FreqHz/2);
                    if notchFreqs(f)==20
                        [bNotch,aNotch] = iirnotch(w0,1/(FreqHz/2));  % 3 dB BW=1 Hz
                    else 
                        [bNotch,aNotch] = iirnotch(w0,4/(FreqHz/2));  % 3 dB BW=4 Hz
                    end
                    
                    rawtrial=filter(bNotch,aNotch,rawtrial,[],1);
                end
            end
            

            if nanArtifacts
                
                % find which element of CndTiming struct matches this
                % condition
                CndTimingIndx=[];
                for ci=1:numel(CndTiming);
                    if CndTiming(ci).cndNmb==TimeLine(t).cndNmb
                        CndTimingIndx=ci;
                    end
                end
                if isempty(CndTimingIndx)
                    error('Can''t find CndTiming variable matching this condition');
                end
                
                thisCondIndx=CndTimingIndx;
                thisPreludeDurSamples=round(CndTiming(thisCondIndx).preludeDurSec*FreqHz);
                thisPostludeDurSamples=round(CndTiming(thisCondIndx).postludeDurSec*FreqHz);
                thisTrialSamples=round(CndTiming(thisCondIndx).nmbTrialSteps*CndTiming(thisCondIndx).stepDurSec*FreqHz);
                if thisPreludeDurSamples+thisPostludeDurSamples+thisTrialSamples~=size(rawtrial,1)
                    warning('unable to verify correctness of artifact rejection');
                end
                
                nEpochsFound=CndTiming(thisCondIndx).nmbTrialSteps;
                if thisPreludeDurSamples>0, nEpochsFound=nEpochsFound+1; end;
                if thisPostludeDurSamples>0, nEpochsFound=nEpochsFound+1; end;
       
                if nEpochsFound==size(IsEpochOK,1)
                    
                    epochDuration=round(size(rawtrial,1)/nEpochsFound);
                    
                    [rowE,colE]=find(IsEpochOK==0);
                    for e=1:length(rowE)
                        rowIdx = ((rowE(e)-1)*epochDuration+1):((rowE(e)-1)*epochDuration+epochDuration);
                        rawtrial(rowIdx,colE(e))=NaN;
                    end
                    
                else
                    warning('Unable to perform artifact rejection due to inconsistency between IsEpochOK and rawdata');
                end
                
            end
            %% remove prelude/postlude
            rawtrial=rawtrial(thisPreludeDurSamples+1:end,:);
            rawtrial=rawtrial(1:end-thisPostludeDurSamples,:);
            
            %% baseline correction
            rawtrial=rawtrial-repmat( nanmean(rawtrial(bslIndx,:),1) , size(rawtrial,1) ,  1 )  ;
            
            %% convert to uV
            rawtrial=rawtrial*10^6;
            
            %% fill in cell array
            rawdata{TimeLine(t).cndNmb}=cat(3,rawdata{TimeLine(t).cndNmb},rawtrial);
            
        end
    end
end

cellData=rawdata;
cd(curDir);
end


%% 
function filelist = subfiles(inputName,incl_path)
    if nargin < 1
        templist = dir;
    else
        templist = dir(inputName);
    end
    if nargin < 2
        incl_path = false;
    else
    end
    tempIdx = strfind(inputName,'/');
    if isempty(tempIdx)
        curDir = pwd;
    else
        curDir = inputName(1:max(tempIdx));
    end
    num_folders = 0;
    for t=1:length(templist)
        if templist(t).isdir ==0 && ~strcmp(templist(t).name,'.') && ~strcmp(templist(t).name,'..')  
            num_folders = num_folders+1;
            if incl_path
                filelist(num_folders,:) = {[curDir,templist(t).name]};
            else
                filelist(num_folders,:) = {templist(t).name};
            end
        else
        end
    end
    if num_folders == 0
        filelist = false;
    else
    end
end
