function[Out,SSX,SSfilt]=BPKF_BandPower_Ratio(X,bandSet,Fs,nWindow)
%% Remark: recommended to drop tail ends later
%% Optionally use [0,f] or [-inf,f] for lowpass, [f,inf] for highpass

X=zscore(X')';
SSX=convn(X.^2,ones(1,nWindow)/nWindow,'same');
if isinf(bandSet(2))
    if nargout<3
    Out=1-(convn((X-highpass(X',bandSet(1),Fs)').^2,ones(1,nWindow)/nWindow,'same')./SSX);
    else
    SSfilt=convn((X-highpass(X',bandSet(1),Fs)').^2,ones(1,nWindow)/nWindow,'same');
    Out=1-(SSfilt./SSX);
    end
elseif isinf(bandSet(1))||bandSet(1)==0
    if nargout<3
    Out=1-(convn((X-lowpass(X',bandSet(2),Fs)').^2,ones(1,nWindow)/nWindow,'same')./SSX);
    else
    SSfilt=convn((X-lowpass(X',bandSet(2),Fs)').^2,ones(1,nWindow)/nWindow,'same');
    Out=1-(SSfilt./SSX);
    end
else
    if nargout<3
    Out=1-(convn((X-bandpass(X',bandSet,Fs)').^2,ones(1,nWindow)/nWindow,'same')./SSX);
    else
    SSfilt=convn((X-bandpass(X',bandSet,Fs)').^2,ones(1,nWindow)/nWindow,'same');
    Out=1-(SSfilt./SSX);
    end
end