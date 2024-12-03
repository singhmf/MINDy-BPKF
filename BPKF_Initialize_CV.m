function[ParStr,BatchSchedule]=BPKF_Initialize_CV(ParStr)


findEmpty=@(xx,yy)(~isfield(xx,yy)||isempty(xx.(yy)));
if findEmpty(ParStr,'CV')
    ParStr.CV.mark=[];
    ParStr.CV.rate=0;
else
    if findEmpty(ParStr.CV,'Size')
        disp('Using BatchSz as unspecified CV size')
        ParStr.CV.Size=ParStr.BatchSz;
    end
end

BatchSchedule=1:ParStr.NBatch;

if ParStr.CV.rate~=0
    for jj=(ParStr.NBatch:-1:1)
        if mod(jj,ParStr.CV.rate)==0
            BatchSchedule=[BatchSchedule(1:jj) -jj/ParStr.CV.rate BatchSchedule((1+jj):end)];
        end
    end
end

end


