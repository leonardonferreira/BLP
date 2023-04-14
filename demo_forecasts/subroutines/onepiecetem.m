function [x_a,fx]=onepiecetem(mean,sd,dof,st_v,grid,n)

%
% Format: {fx}=onepiecetem(mean,sd,dof,st_v,grid,n);
% Inputs:  mean: the mean of distribution
%          sd:   the standard deviation of the distribution
%          dof:  number of the degrees of freedom
%          st_v: starting value of the x axis
%          gird: the increment by which the x axis is increased
%          n:    the end point of the x axis; this would be the row dimension of the output
% Output:  fx:   an nx2 matrix. The first column is the grid of the x axis, and the second
%                column contains the t-distribution values, normalised to sum to 100
%



var=sd^2; x_a=(st_v:grid:n)'; 
fx=zeros(size(x_a,1),1);
for i=1:size(fx,1); 
    fx(i,1)=mstudent(x_a(i,1),mean,var,dof); 
end;







