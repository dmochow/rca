function Y = rcaProject(data,W)
% Y=RCAPROJECT(DATA3D,W)
%
% data3D: array of epoched data (samples x channels x trials)
% W: matrix of RCA projection vectors (channels x component)
%
% Y: RCA projected data (samples x components x trials)
%
% throughout data, NaNs indicate missing data values (rejected artifacts)
%
% (c) Jacek P. Dmochowski, 2014


if ~iscell(data)
    [nSamples]=size(data,1);
    [nElectrodes]=size(data,2);
    [nTrials]=size(data,3);
    if nSamples<nElectrodes, warning('Number of samples must be greater than or equal to the number of electrodes'); end;
    if nElectrodes~=size(W,1), error('dimension 2 of data must match dimension 1 of W'); end
    nComp=size(W,2);
    
    Y=zeros(nSamples,nComp,nTrials);
    for comp=1:nComp
        Y(:,comp,:)= squeeze ( nansum ( data .* repmat(W(:,comp).',[nSamples 1 nTrials]) , 2 ) )  ;
    end
    
else  % cell mode
    [nCond,nSubjects]=size(data);
    nComp=size(W,2);
    Y=cell(nCond,nSubjects);
    for c=1:nCond
        for s=1:nSubjects

            data3D=data{c,s};
            [nSamples,nElectrodes,nTrials]=size(data3D);
            if nSamples<nElectrodes, warning('Number of samples must be greater than or equal to the number of electrodes'); end;
            if nElectrodes~=size(W,1), error('dimension 2 of data must match dimension 1 of W'); end
            
            for comp=1:nComp
                Y{c,s}(:,comp,:)= squeeze ( nansum ( data3D .* repmat(W(:,comp).',[nSamples 1 nTrials]) , 2 ) )  ;
            end
                                   
        end
    end
    
end
