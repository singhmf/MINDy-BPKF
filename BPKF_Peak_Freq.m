function[fOut,pOut,boxPSD,Fmark]=BPKF_Peak_Freq(x,Fband,Fs,fRes,nSmth)
x=x-mean(x,2);

N=size(x,2);
szX=size(x);
sz2=szX;sz2(2)=floor(N/2)+1;
psdx=zeros(sz2);
xdft = fft(x,[],2);
psdx(:,:,:) = (1/(Fs*N)) * abs(xdft(:,1:N/2+1,:)).^2;
%psdx(2:end-1) = 2*psdx(2:end-1);

freq = 0:(Fs/N):Fs/2;

Fmark=Fband(1):fRes:Fband(2);

boxPSD=zeros(size(psdx,1),numel(Fmark));
if nSmth~=1
freq=convn(freq,ones(1,nSmth)/nSmth,'valid');
boxPSD=convn(boxPSD,ones(1,nSmth)/nSmth,'valid');
end
for ii=1:numel(Fmark)
    %% Discretization on-center
    goodFreq=abs(freq-Fmark(ii))<=(fRes/2);
    boxPSD(:,ii)=mean(psdx(:,goodFreq),2);
end
[pOut,t2]=max(boxPSD,[],2);
fOut=Fmark(t2);
end

