function[mOut,bOut,predOut,residOut]=BPKF_TheilSen_Regress(Xvec,Yvec)

mOut=median(OffDiag((Yvec-Yvec'))./OffDiag((Xvec-Xvec')),'all');
bOut=median(Yvec-mOut*Xvec);
if nargout>2
    predOut=mOut*Xvec+bOut;
end
if nargout>3
    residOut=Yvec-predOut;
end
end