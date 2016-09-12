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


function [X,q]=fPRm(A,alpha,M,Ngrid,I)

% The following program takes in input a cell A of 2 sparse adjacency 
% matrices (ordered layers of the duplex networks A{1},A{2}) and computes 
%
%  - the Functional Multiplex PageRank X(:,I(i)) of node I(i) for all values of the 
%     influences z=[z_1,z_2,..z_M] of multilinks with multiplicity of
%     overlap respectivelly 1,2,3,...M
%     taking values z=[1,q,q^2,...q^(M-1)] 
%     with q taking  positive values.
%   -alpha is the damping factor (usually taken 0.85)
%
%    -Ngrid determines the number of points in the grid (q_r)
%           with q_r=exp((r-Ngrid*0.5)*0.2) with r=1,2...,Ngrid
%        
%
%    -I indicates the label of the nodes that we want to consider.
%       the default setting is 1:N 
%       but to save memory if is possible to choose a vector I including
%       only a subset of the total labels.
%  Usage:
% to represent the Functional Multiplex Pagerank of node I(i) plot
% contour(q,X(:,I(i)))

% Notes:
% The adjacency matrices A{1} and A{2},.. A{M} must have the same dimension and are
% unweighted
% For the  adjacency matrices we adopt the following convention
% A{p}(i,j)=1 if and only if j points to in layer p otherwise A{p}(i,j)=0




for i1=1:Ngrid,
    q(i1)=exp((i1-Ngrid*0.5)*0.2);
    for i=1:M,
        z(i)=(q(i1))^(i-1);
    end
    
[x]=functionalPageRank_multiplicity(A,alpha,z,M);
for i=1:numel(I),
X(i1,i)=x(I(i));
end

 end

end