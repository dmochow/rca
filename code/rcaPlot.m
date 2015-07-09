function [h,hlg,grpData,grpDataSem] = rcaPlot( rcaData , A , grpIndx , t , colors, labels , bslSample , locfile)
%[h,hlg,grpData,grpDataSem] = rcaPlot( rcaData , A , grpIndx , t , colors, labels , bslSample )
%   show forward models of components and corresponding time courses
%   THIS FUNCTION STILL UNDER DEVELOPMENT: OMEGA VERSION

if nargin<8, locfile='GSN-HydroCel-128.sfp'; end
if nargin<7, bslOn=0; else bslOn=1; end;

if iscell(rcaData)
    nCond=size(rcaData,1);
    nGroups=length(grpIndx);
    [nSamples,nComp,~]=size(rcaData{1,1});
    if max([grpIndx{:}])>nCond, error('JD: grpIndx element exceeds number of conditions'); end
else
    %[nSamples,nComp,nTrials]=size(rcaData);
    error('JD: so far, only support for cell data input');  
end

grpData=zeros(nSamples,nComp,nGroups);
grpDataSem=zeros(nSamples,nComp,nGroups);
for g=1:nGroups
    tmp=cat(3,rcaData{grpIndx{g},:});
    
    if bslOn
        bsl=tmp(bslSample,:,:);
        tmpBsl=tmp-repmat(bsl,[size(tmp,1) 1 1]);
        grpData(:,:,g)=nanmean(tmpBsl,3);  % baseline
    else
        grpData(:,:,g)=nanmean(tmp,3);  % no baseline
    end
    
    grpDataSem(:,:,g)=nanstd(tmp,[],3)/sqrt(size(tmp,3));  
end



figure
wd=0.2; 
wd_=0.7;
h=zeros(nComp*2,1);
hlg=zeros(nComp,1);
clr=colors; 
for c=1:nComp
    
    h((c-1)*2+1)=subplot(nComp,2,(c-1)*2+1); topoplot(A(:,c), locfile,'numcontour',0); axis off; 
    pos=get(h((c-1)*2+1),'Position');
    set(h((c-1)*2+1),'Position',[0 pos(2) wd pos(4)])
    
    h(c*2)=subplot(nComp,2,c*2); 
    %plot(t,squeeze(grpData(:,c,:))); 
    for g=1:nGroups
       hsh(g)=shadedErrorBar(t,squeeze(grpData(:,c,g)),squeeze(grpDataSem(:,c,g)),clr(g),1); hold on;
    end
    pos=get(h(c*2),'Position');
    set(h(c*2),'Position',[wd pos(2) wd_ pos(4)]);
    xlim([t(1) t(end)])
    %if c==nComp
    if 1
        %hlg=legend(labels);
        hlg(c)=legend([hsh(:).patch],labels);
        %lgPos=get(hlg,'Position');
        %set(hlg,'Position',[lgPos(1)+0 lgPos(2) lgPos(3) lgPos(4)]);
        set(hlg(c),'Box','off');
        set(hlg(c),'Location','NorthWest','orientation','horizontal');
    end;
end


end

