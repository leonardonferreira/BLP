function [x_a,fx]=onepiecenem(mean,sd,st_v,grid,n)

% A procedure to graph a normal distribution. 
% 
% Format:  {fx}=onepiecenem(mean,sd,st_v,grid,n)
% Inputs:  mean: the mean of the distribution
%          sd:   standard deviation of the distribution
%          st_v: the starting value of the x-axis
%          grid: the grid of the x-axis, by how much it increases in every counter
%          n:    the top end of the x-axis, the number of rows of the distribution
%
% Outputs: nx2 matrix. The first column contains the labelling of the x-axis, with the
%          second column having the densities.         
%




x_a=(st_v:grid:n)';
fx=zeros(size(x_a,1),1);
for i=1:size(fx,1); 
    fx(i,1)=normpdf(x_a(i,1),mean,sd);
    %fx(i,1)=pdf('norm',(((x_a(i,1))-mean)/sd),0,1)/sd; 
end;

