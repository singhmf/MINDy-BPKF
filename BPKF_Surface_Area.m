function[PosArea,TriArea]=BPKF_Surface_Area(SourceMesh)
%% Must have fields pos and tri

VecSet=cell(1,3);
for ii=1:3
VecSet{ii}=SourceMesh.pos(SourceMesh.tri(:,ii),:);
end

TriArea=.5*sqrt(sum(cross(VecSet{1}-VecSet{2},VecSet{3}-VecSet{2},2).^2,2));
PosArea=zeros(size(SourceMesh.pos,1),1);

for ii=1:size(SourceMesh.tri,1)
    PosArea(SourceMesh.tri(ii,:))=PosArea(SourceMesh.tri(ii,:))+TriArea(ii);
end


end