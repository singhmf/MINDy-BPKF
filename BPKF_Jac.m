function[Out,dSig]=BPKF_Jac(ooP,Xeval)
%% Evaluate Jacobian
dSig=ooP.Param{4}.*(1-tanh(ooP.Param{4}.*Xeval+ooP.Param{5}).^2);
if size(Xeval,2)~=1
    disp('Evaluating average Jacobian')
    Out=ooP.Param{1}.*mean(dSig,2)'+diag(ooP.Param{2});
else
Out=ooP.Param{1}.*dSig'+diag(ooP.Param{2});
end