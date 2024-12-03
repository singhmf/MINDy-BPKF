function[X1,X2]=BPKF_Split_Cell(Xset)
nX=cellfun(@(xx)(size(xx,2)),Xset);

mm=find(cumsum(nX)==(sum(nX)/2));

if ~isempty(mm)
    X1=Xset(1:mm);
    X2=Xset((mm+1):end);
else
midMark=ceil(sum(nX)/2);
CombSet=[0 cumsum(nX)];
meanMark=find(and(CombSet(1:end-1)<=midMark,...
                  midMark<=CombSet(2:end)));
SplitPoint=midMark-CombSet(meanMark);
X1=[Xset(1:meanMark-1), {Xset{meanMark}(:,1:SplitPoint)}];
X2=[{Xset{meanMark}(:,(SplitPoint+1):end)},Xset(meanMark+1:end)];
end
end