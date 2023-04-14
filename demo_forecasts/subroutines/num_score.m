function [oo]=num_score(density,outturn,st_v,grid)
% 
% st_v - starting value on the left 
% grid - grid size 
% 
% numerical calculation of score 
v=outturn;
nn=indexcat(v,st_v,grid);
oo=density(nn);
if oo==0; oo=0.00001;end;



function [c]=indexcat(v,st_v,grid)
c=0;
i=st_v;
while i<v;
c=c+1;
i=i+grid;
end;


