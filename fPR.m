%++++++++ Functional Multiplex PageRank ++++++++++++++++++++++++++++++++++++++++++++
% 
% This code can be redistributed and/or modified
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%  
% This program is distributed ny the authors in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%  
% If you use this code please cite 
%
% [1]  J. Iacovacci, C. Rahmede, A. Arenas and G. Bianconi, "Functional Multiplex PageRank." 
%        arxiv:1608.06328 (2016)  
%
% (c) Jacopo Iacovacci (mriacovacci@hotmail.it) 
%     Christoph Rahmede (c.rahmede@kit.edu)
%     Ginestra Bianconi (g.bianconi@qmul.ac.uk)  
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [X,theta,phi]=fPR(A,alpha,Ngrid,I)
% The following program takes in input a cell A of 2 sparse adjacency 
% matrices (ordered layers of the duplex networks A{1},A{2}) and computes 
%
%  - the Functional Multiplex PageRank X(:,:,I(i)) of node I(i) for all values of the 
%     influences z=[z^(1,0),z^(0,1),z^(1,1)]
%     taking values z=[sin(theta)cos(phi),sint(theta)sin(phi),Cos(theta)] 
%     with 0<=theta<=pi/0.5 and 0<=phi<=pi/0.5
%   
%   -alpha is the damping factor (usually taken 0.85)
%
%    -Ngrid determines the number of points in the grid (theta_r,phi_s)
%           with theta_r=(pi*0.5/Ngrid)*r with r=0,1,2...,Ngrid
%            and phi_rs(pi*0.5/Ngrid)*s with s=0,1,2...,Ngrid
%
%    -I indicates the label of the nodes that we want to consider.
%       the default setting is 1:N 
%       but to save memory if is possible to choose a vector I including
%       only a subset of the total labels.
%  Usage:
% to represent the Functional Multiplex Pagerank of node I(i) plot
% contour(theta,phi,X(:,:,I(i)))

% Notes:
% The adjacency matrices A{1} and A{2} must have the same dimension and are
% unweighted
% For the  adjacency matrices we adopt the following convention
% A{1}(i,j)=1 if and only if j points to in layer 1 otherwise A{1}(i,j)=0
% A{2}(i,j)=1 if and only if j points to in layer 2 otherwise A{2}(i,j)=0



for i1=1:(Ngrid+1),
    theta(i1)=((pi*0.5)/Ngrid)*(i1-1);
    for i2=1:(Ngrid+1),
        phi(i2)=((pi*0.5)/Ngrid)*(i2-1);
        z(1)=sin(theta(i1))*cos(phi(i2));
        z(2)=sin(theta(i1))*sin(phi(i2));
        z(3)=cos(theta(i1));

	[x]=functionalPageRank_duplex(A,alpha,z);
	for i=1:numel(I),
	X(i1,i2,i)=x(I(i));
    end
    end
end