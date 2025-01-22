function[Out]=OffDiag(X)
X=squeeze(X);
nX=size(X);
if nX(1)~=nX(2)
    error('X should be square')
end
nX0=nX;nX0([1 2])=1;
if ~ismatrix(X)
%    disp('Warning X is not a matrix, taking off for the first 2 dims only');
Out=X(repmat(eye(nX(1)),nX0)==0);
else
Out=X(eye(nX)==0);
end
end