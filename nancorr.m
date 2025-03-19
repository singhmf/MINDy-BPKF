function[Out,OutP]=nancorr(X,Y,varargin)
%% Uses all rows
BadInd=or(isinf(X+Y),isnan(X+Y));
X(BadInd)=[];Y(BadInd)=[];
if ~isempty(X)
[Out,OutP]=corr(X,Y,'Rows','all',varargin{:});
else
    Out=nan;
    OutP=nan;
end
end