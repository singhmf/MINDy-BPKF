function[]=BPKF_Plot_Err_Huber(Out,SmthSpace,doCV,varargin)

if nargin==1
    disp('Assuming no smoothing kernel')
    SmthSpace=1;
end

if nargin==2
if isfield(Out,'kalCV')
    if isfield(Out,'batchCV')
        batchCV=Out.batchCV;
    else
        disp('Assuming CV rate is 500')
        batchCV=(1:size(Out.kalCV_Huber,2))*500;
    end
    doCV=true;
else
    doCV=false;
end
else
    doCV=strcmpi(doCV(1),'y');
end


if numel(varargin)==1
    varargin=repelem(varargin,1,2);
end

figure
subplot(1,2,1)
plot(SmthSpace:size(Out.kalE_Huber,2),convn(Out.kalE_Huber,ones(1,SmthSpace)/SmthSpace,'valid')');
set(gca,'ColorOrderIndex',1)
if doCV
    hold on;plot(Out.batchCV,Out.kalCV_Huber,'.')
end
if ~isempty(varargin)&& strcmpi(varargin{1},'y')
    set(gca,'yscale','log')
end
title('Kal Error')
pbaspect([1 1 1])

subplot(1,2,2)
plot(SmthSpace:size(Out.recE_Huber,2),convn(Out.recE_Huber,ones(1,SmthSpace)/SmthSpace,'valid')')
set(gca,'ColorOrderIndex',1)
if doCV
    hold on;plot(Out.batchCV,Out.recCV_Huber,'.')
end
if ~isempty(varargin)&& strcmpi(varargin{2},'y')
    set(gca,'yscale','log')
end
title('Rec Error')
pbaspect([1 1 1])
end