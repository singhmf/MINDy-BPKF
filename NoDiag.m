function[Out]=NoDiag(M,varargin)
%% Varargin allows a block size
%% If greater than 2D always performs along the [1 2] axis
M=squeeze(M);
if isempty(varargin)
    nBlock=1;
elseif ceil(varargin{1})==varargin{1}&&varargin{1}>0
    nBlock=varargin{1};
else
    error('Block size should be empty or natural')
end
nM=size(M);
Out=M.*(1-kron(eye(nM([1 2])/nBlock),ones(nBlock)));
end