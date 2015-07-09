function [Wsub,A,W,dGenSort,K] = rcaTrain(Rxx,Ryy,Rxy,K,C)
% [WSUB,A,W]=RCATRAIN(RXX,RYY,RXY,[K],[C])
%
% compute spatial filters maximizing reliability across trials given precomputed auto- and cross-covariance matrices 
%
% Rxx: covariance matrix of "data record 1" (see Dmochowski et al., Front Hum Neurosci 2012)
% Ryy: covariance matrix of "data record 2"
% Rxy: cross-covariance matrix of data records 1 and 2
% K: number of autocovariance dimensions to diagonalize (defaults to the number of eigenvalues explaining 60% of total power)
% C: number of components to return (defaults to 3)
%
% Wsub: channel x component matrix of top C reliability-maximizing projection [Wsub(:,1),...,Wsub(:,C)]
% A: channel x component matrix of corresponding forward models [A(:,1),...,A(:,C)]
% W: channel x component matrix of all RCA projections [W(:,1),...,W(:,nChannels)]
% dGenSort: vector of sorted (ascending) generalized eigenvalues (value conveys correlation coefficient)
%
% (c) Jacek P. Dmochowski, 2014
 
% temporary: ask user what to do here
% nChannels=size(Rxx,1);
% if sum(isnan(Rxx(:)))~=0
%     colsum=sum(isnan(Rxx),1);
%     badIndxXX=find(colsum==size(Rxx,1));
% end
% 
% if sum(isnan(Ryy(:)))~=0
%     colsum=sum(isnan(Ryy),1);
%     badIndxYY=find(colsum==size(Ryy,1));
% end
% 
% if sum(isnan(Rxy(:)))~=0
%     colsum=sum(isnan(Rxy),1);
%     badIndxXY=find(colsum==size(Rxy,1));
% end
% 
% goodIndx=setdiff(1:size(Rxx,1),unique([badIndxXX badIndxYY badIndxXY]));
% Rxx=Rxx(goodIndx,goodIndx);
% Ryy=Ryy(goodIndx,goodIndx);
% Rxy=Rxy(goodIndx,goodIndx);
if sum(isnan(Rxx(:)))~=0 || sum(isnan(Ryy(:)))~=0 || sum(isnan(Rxy(:)))~=0
    warning('Covariance matrices contain NaNs: setting to zero...');
    Rxx(isnan(Rxx))=0; Ryy(isnan(Ryy))=0; Rxy(isnan(Rxy))=0; 
end
%Rxx(isnan(Rxx))=0; Ryy(isnan(Ryy))=0; Rxy(isnan(Rxy))=0; 
% end temporary

if nargin<5 || isempty(C), C=3; end
if nargin<4 || isempty(K)
    if nargin<3, error('At least 3 arguments required'); end
    if sum(size(Rxx)~=size(Ryy)), error('Rxx and Ryy must have the same size'); end
    if sum(size(Rxx)~=size(Rxy)), error('Rxx and Rxy must have the same size'); end
    [~,eigDiag]=eig(Rxx+Ryy); %eigvals=sort(diag(eigDiag),'descend'); propExpl=cumsum(eigvals)/sum(eigvals);
    [~,indxKnee] = knee_pt(diag(eigDiag),1:length(diag(eigDiag)));
    K=length(diag(eigDiag))-indxKnee;
    %thresh=0.6; % conservative for now
    %K=find(propExpl>thresh,1,'first');
    fprintf('Using %d bases to diagonalize pooled autocovariance \n',K+1);
end

[Vpool,Dpool]=eig(Rxx+Ryy); 
[dPool,sortIndx]=sort(diag(Dpool),'ascend');
Vpool=Vpool(:,sortIndx);  % in case eigenvectors/eigenvalues not sorted
dPool=dPool(end-K:end);
Rw=Vpool(:,end-K:end)*diag(1./dPool)*Vpool(:,end-K:end)'*(Rxy+Rxy');  % regularized inverse

[Vgen,Dgen]=eig(Rw);  % compute generalized eigenvalues/eigenvectors
dGen=diag(Dgen);
%[dGenSort,b]=sort(real(dGen)); 
[dGenSort,b]=sort(abs(dGen));  % 11/18/14: complex RCA (?)
W=Vgen(:,b(end:-1:1)); % in case not sorted
%W=real(W);  % ignore small imaginary component

Wsub=W(:,1:C);  % return only selected number of components
Rpool=0.5*(Rxx+Ryy);  
A=Rpool*Wsub*inv(Wsub'*Rpool*Wsub);  % return forward models (see Parra et al., 2005, Neuroimage)

% % put back to full size
% 
% WW=zeros(nChannels,nChannels);
% WW(goodIndx,goodIndx)=W;
% 
% AA=zeros(nChannels,C);
% AA(goodIndx,:)=A;
% 
% Wsub_=zeros(nChannels,C);
% Wsub_(goodIndx,:)=Wsub;
% 
% W=WW; A=AA; Wsub=Wsub_;

