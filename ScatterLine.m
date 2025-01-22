function[Out,r]=ScatterLine(X,Y,varargin)
ttTemp=and(and(~isnan(X(:)),~isnan(Y(:))),and(~isinf(X(:)),~isinf(Y(:))));
Out=scatter(X(ttTemp),Y(ttTemp),varargin{:});
lsline
if isinteger(X)
    X=single(X);
end
if isinteger(Y)
    Y=single(Y);
end
r=(nancorr(X(:),Y(:)));
title(strcat('r=',num2str(r)));
grid on
end