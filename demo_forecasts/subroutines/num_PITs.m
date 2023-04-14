function [PIT]=num_PITs(density,outturnv,st_v,grid,n)

% Format: {PIT}=num_PITs(density,outturnv,st_v,grid,n);
% Inputs: density:  an n times 1 vector of density
%         outturnv: the realisation of the observation we are trying to forecast
%         st_v:     the starting value of the x-axis
%         grid:     the increment by which the x-axis increases
%         n:        the dimension of the x-axis
% Output: score: a scaler of PIT
%
%

   xx=(st_v:grid:n)'; % x-axis

   % numerical calculation of pit's 
   PIT=num_pit(xx,density,outturnv,st_v,grid); 
   

function [ zit ] = num_pit(xx,comb,outturn,st_v,grid)
%  numerical calculation of PITS

v1=outturn-(grid/10);
v2=outturn+(grid/10);
v=[v1 ;v2 ];
v=outturn;
zitt=comb.*grid;
nn=indexcat(xx,v,st_v,grid);
zit=sum(zitt(1:nn));

function [c]=indexcat(xx,v,st_v,grid)
c=0;
i=st_v;
while i<v;
c=c+1;
i=i+grid;
end;












