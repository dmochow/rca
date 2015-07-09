function dataOut= baselineVolume(dataIn,bslSample )
%dataOut= baselineVolume(dataIn,bslSample )
%   shift each trial (dim 3) to a consistent baseline level
if ndims(dataIn)~=3, error('JD: data must have three dimensions'); end
if size(dataIn,1)<size(dataIn,2), dataIn=permute(dataIn,[2 1 3]); warning('JD: permuting dimensions 1 and 2'); end
if bslSample<1 || bslSample>size(dataIn,1), error('JD: out of range baseline sample'); end

bsl=dataIn(bslSample,:,:);
dataOut=dataIn-repmat(bsl,[size(dataIn,1) 1 1]);

end

