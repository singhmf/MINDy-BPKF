function[minDist]=BPKF_Mesh_Distances(SourceMesh,doReScale,testChan)
%% SourceMesh contains .pos,.tri
%% qExpC (cell or array) contains exponents for kernels
%% doReScale (y/n) denotes whether to rescale distances so max range=1 in largest directions or use native units.
%% arOrder=1: exponential    arOrder=2: gaussian



mm.vertices=SourceMesh.pos;

if strcmpi(doReScale(1),'y')
    mm.vertices=mm.vertices/max(range(mm.vertices,1));
end

mm.faces=SourceMesh.tri;
nVert=size(mm.vertices,1);

%% Get 2norm distances
%DistMat=zeros(nVert);
%for ii=1:3
%    DistMat=DistMat+(mm.vertices(:,ii)-mm.vertices(:,ii)').^2;
%end
%DistMat=sqrt(DistMat);

%% Make surface-connection graph
%DistMat=sparse(nVert);

compMat=[1 -1 0; 1 0 -1; 0 1 -1];
zz=false(nVert);
for ii=1:size(mm.faces,1)
    pp=mm.faces(ii,:);
    zz(pp([2 3]),pp(1))=true;
    zz(pp(2),pp(3))=true;
    zz(pp(3),pp(2))=true;
end
DistMat=sparse((zz));
for ii=1:size(mm.faces,1)
    tmp=sum((compMat*mm.vertices(pp,:)).^2,2);
    DistMat(pp(1),pp([2 3]))=tmp([1 2]); %#ok<SPRIX> 
    DistMat(pp(2),pp(3))=tmp(3); %#ok<SPRIX> 
end
DistMat=DistMat+DistMat';
%% Minimum distance along surface

if nargin<3
    minDist=distances(graph(DistMat));
else
    minDist=distances(graph(DistMat),testChan);
end

end
