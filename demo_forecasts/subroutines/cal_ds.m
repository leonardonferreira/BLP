function calds=cal_ds(fyds,fmds,nds)
% -- Construct Calendar Sequences for Plotting, etc -- @
%    fyds=1987; fmds=4; 

calds=zeros(nds,2); calds(1,1)=fyds; calds(1,2)=fmds;
yr=fyds; mt=fmds;
for i=2:nds;
    mt=mt+1;
    if mt > 4; mt=1; yr=yr+1; end; 
    calds(i,1)=yr; calds(i,2)=mt;
end; 