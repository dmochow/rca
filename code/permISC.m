function [ISC, crit_ISC, pval, est_alpha] = permISC(data,W,Pinfo)
%[ISC, crit_ISC, pval, est_alpha] = permISC(data,W,Pinfo)
%
% Calculates intra subject correlations (ISC) on data projected using given
% spatial filters, W. A syntaxt has been used to complement the RCA scrips
% supplied by Jacek Dmochowski. The script is inspired by papers from the
% same author (2012 and 2014).
%
% (c) Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, Cognitive systems - November 2014
%
% Inputs
% data  - A 3-D array (samples x channels x trials) % data can either be a
%         3-D array (samples x channels x trials), or a cell array
%         (conditions x subjects) of 3-D data arrays (as above).
% W     - Channel x component matrix of all RCA projections [W(:,1),...
%         ,W(:,nChannels)].
% Pinfo - [optional. Can also be set as "[]" for default settings]. Struct
%         containing fields that define properties for correlation and
%         permutation test. It is possible to only define some of the
%         fields.
%   .win - Number of samples in each window to be correlated.
%  .step - Number of samples the window is moved between correlations.
% .Nperm - Number of permutations conducted for a test of significance.
% .alpha - Desired alpha level for permutation test.
% .tail  - 1 for upper one-tailed test. 2 for two-tailed test.
%
% Outputs
% ISC       - Array containing ISC calculated for each window and component.
% crit_ISC  - Estimated upper and lower critical correlation coefficients
%             for given alpha level.
% pval      - Estimated p-values.
% est_alpha - If Nperm is low, it might not be possible to use the precise
%             alpha value. Instead a value as close as possible is used.
%             This value is estimated in est_alpha.

%% Checking inputs
if nargin<3
    Pinfo = [];
end
if ~isfield(Pinfo,'step'),	step = 128; else step = Pinfo.step; end
if ~isfield(Pinfo,'win'),	win = 5 * step; else win = Pinfo.win; end
if ~isfield(Pinfo,'Nperm'),	Nperm = 0; else Nperm = Pinfo.Nperm; end
if ~isfield(Pinfo,'alpha'),	alpha = 0.01; else alpha = Pinfo.alpha; end
if ~isfield(Pinfo,'tail'),	tail = 1; else tail = Pinfo.tail; end

if ~iscell(data),
    Ncond=1; Nsubs=1; subjRange=1; condRange=1; data={data};
else
    Ncond=size(data,1);
    if ~isfield(Pinfo,'condRange'),	condRange=1:Ncond;
    else condRange = Pinfo.condRange; end
    Nsubs=size(data,2);
    if ~isfield(Pinfo,'subjRange'),	subjRange=1:Nsubs;
    else subjRange = Pinfo.subjRange; end
    data=data(condRange,subjRange);
end

Nc = size(W,2);

%% Calculating permutation order
% Assuming different subjects and conditions have the same number of samples
N = size(data{1,1},1);
Nwin = floor((N-win)/step);
ISC = zeros(Nc,Nwin);
if tail == 1
    crit_ISC = zeros(Nc,Nwin);
elseif tail == 2
    crit_ISC = zeros(Nc,Nwin,2);
else
    error(['Wrong value of tail given. Should be 1 or 2, was set to tail='...
        num2str(tail) '.'])
end
pval = zeros(Nc,Nwin);
est_alpha = zeros(Nc,Nwin);

permorder = cell(Ncond,Nsubs);
for cnd=1:Ncond
    for sub=1:Nsubs
        Ntrials = size(data{cnd,sub},3);
        Ncomb = Ntrials*(Ntrials-1)/2;% number of trial combinations
        
        % Preallocating the time shift of each trial pair for each permutation
        permorder{cnd,sub} = zeros(Ncomb,Nperm,'uint16');
        for n = 1:Nperm
            for nn = 1:Ncomb
                permorder{cnd,sub}(nn,n) = round(rand*Nwin-.5);
            end
        end
    end
end
%% Calculating ISC with permutations for one window at a time.
tstart = tic;
for w_ind = 1:Nwin
    [Rxx, Ryy, Rxy] = calcRs(data,w_ind,step,win,Nwin,[],[],0);
    ISC(:,w_ind) = calcISC(Rxx,Ryy,Rxy,W);
    pISC = zeros(Nc,Nperm);
    for p = 1:Nperm
        if rem(p,500)==0
            fprintf('Window no.: %d, permutation no.: %d, Time spent.: %d s.\n',...
                w_ind,p,round(toc(tstart)))
        end
        [Rxx, Ryy, Rxy] = calcRs(data,w_ind,step,win,Nwin,permorder,p,1);
        pISC(:,p) = calcISC(Rxx,Ryy,Rxy,W);
    end
    
    
    if Nperm>0
        % Significance testing
        for c = 1:Nc
            if ISC(c,w_ind)>0 || tail==1
                pval(c,w_ind)=mean(pISC(c,:) > ISC(c,w_ind))*tail;
            else
                pval(c,w_ind)=mean(pISC(c,:) < ISC(c,w_ind))*tail;
            end
            crit_ISC(c,w_ind,1) = prctile(pISC(c,:),100*(1-alpha/tail));
            est_alpha(c,w_ind) = mean(pISC(c,:) >= crit_ISC(c,w_ind,1));
            
            if tail == 2
                crit_ISC(c,w_ind,2) = prctile(pISC(c,:),100*alpha/tail);
                est_alpha(c,w_ind) = est_alpha(c,w_ind) + mean(pISC(c,:) <= crit_ISC(c,w_ind,2));
            end
        end
    end
end

if Nperm == 0
    crit_ISC = [];
    est_alpha = [];
    pval = [];
end

end

function ISC = calcISC(Rxx,Ryy,Rxy,W)
ISC = (sum((W'*Rxy)'.*W) ./ sqrt(sum((W'*Rxx)'.*W).*sum((W'*Ryy)'.*W)))';
end

function [Rxx, Ryy, Rxy] = calcRs(data,w_ind,step,win,Nwin,permorder,p_ind,perm)
% Finding sample index for the given window.
idx = 1+step*(w_ind-1):step*(w_ind-1)+win;
Nsample = length(idx);

Ncond = size(data,1);
Nsubs = size(data,2);
[~,Nelectrode,~] = size(data{1,1});
% Preallocating
sumXX = zeros(Nelectrode,Nelectrode);
nPointsInXX = zeros(Nelectrode,Nelectrode);
sumYY = zeros(Nelectrode,Nelectrode);
nPointsInYY = zeros(Nelectrode,Nelectrode);
sumXY = zeros(Nelectrode,Nelectrode);
nPointsInXY = zeros(Nelectrode,Nelectrode);

for cnd=1:Ncond
    for sub=1:Nsubs
        Ntrials = size(data{cnd,sub},3);
        pindx = combnk(1:Ntrials,2);
        
        thisVolume = data{cnd,sub}(idx,:,:);
        concat1 = permute(thisVolume(:,:,pindx(:,1)),[2 1 3]);
        
        if perm
            % Permuting the second view by timeshift given in permorder
            concat2 = zeros(Nelectrode,Nsample,size(pindx,1));
            for comb = 1:size(pindx,1)
                idx_perm = calcPermidx(w_ind,step,win,Nwin,permorder{cnd,sub}(comb,p_ind));
                trial = pindx(comb,2);
                concat2(:,:,comb) = data{cnd,sub}(idx_perm,:,trial)';
            end
        else
            concat2 = permute(thisVolume(:,:,pindx(:,2)),[2 1 3]);
        end
        
        % Lucas' trick to ensure that Rxx=Ryy
        concatX = [concat1(:,:) concat2(:,:)];
        concatY = [concat2(:,:) concat1(:,:)];
        
        % offset by mean local to this subject-condition
        concatX = concatX-repmat(nanmean(concatX,2),[1 size(concatX,2)]);
        concatY = concatY-repmat(nanmean(concatY,2),[1 size(concatY,2)]);
        
        nPointsInXX = nPointsInXX + double(~isnan(concatX))*double(~isnan(concatX))';
        nPointsInYY = nPointsInYY + double(~isnan(concatY))*double(~isnan(concatY))';
        nPointsInXY = nPointsInXY + double(~isnan(concatX))*double(~isnan(concatY))';
        
        concatX(isnan(concatX)) = 0;
        concatY(isnan(concatY)) = 0;
        
        sumXX = sumXX + concatX*concatX';
        sumYY = sumYY + concatY*concatY';
        sumXY = sumXY + concatX*concatY';
    end
end
Rxx = sumXX ./ nPointsInXX;
Ryy = sumYY ./ nPointsInYY;
Rxy = sumXY ./ nPointsInXY;
end

function idx_perm = calcPermidx(w_ind,step,win,Nwin,permorder)
% Finding sample index for permuted window. The timeshift is given in
% permorder.
offset = w_ind + permorder - 1;
w_perm = rem(offset,Nwin) + 1;
idx_perm = 1+step*(w_perm-1):step*(w_perm-1)+win;
end
