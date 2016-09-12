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
%        arXiv:1608.06328 (2016).  
%
% (c) Jacopo Iacovacci (mriacovacci@hotmail.it) 
%     Christoph Rahmede (c.rahmede@kit.edu)
%     Ginestra Bianconi (g.bianconi@qmul.ac.uk)  
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [x]=functionalPageRank_multiplicity(A,alpha,z,M)

% The following program takes in input a cell A of M sparse adjacency 
% matrices (ordered layers of the duplex networks A{1},A{2},...A{M}) and computes 
%
%  - the Functional Multiplex PageRank x_i of node i given the vector of
%     influences z=[z_1,z_2,z_3...z_M] of multilinks with multiplicity of
%     overlap respectivelly 1,2,...M.
%  
%   -alpha is the damping factor (usually taken 0.85)
%
%   -M is the total number of layers of the multiplex network
%
% Options:
%
% The influences z=[z_1,z_2,z_3...z_M] of multilinks with multiplicity of
%     overlap respectivelly 1,2,...M.
% Notes:
% The adjacency matrices A{1},A{2}... A{M} must have the same dimension and are
% unweighted
% For the  adjacency matrices we adopt the following convention
% A{p}(i,j)=1 if and only if j points to in layer p otherwise A{p}(i,j)=0



v_quadratic_error=0.0001;

Am=zeros(size(A{1}));
for n=1:M,
    Am=Am+A{n};
end
for n=1:M,
    B{n}=(Am==n);
end


x0=zeros(1,max(size(B{1}))); 

for jj=1:M
       
x0=x0+z(jj)*((sum(B{jj},1)>0)+(sum(B{jj},2)>0)');
end
x0=(x0>0);
norma=nnz(x0);
x0=x0/norma;
x0=x0';
D=zeros(1,max(size(B{1})));    
for jj=1:M,
D=D+sum(B{jj},1)*z(jj);
end
D=D+(D==0);
D=ones(1,max(size(B{1})))./D;
n=numel(D);
D=spdiags(D(:),0,n,n);
l=zeros(size(sum(B{1},1)));
for n=1:M,
    l=l+z(n)*(sum(B{n},1));
end
l=l>0;
jump=alpha*l';

Am2=zeros(size(B{1}));
for n=1:M,
    Am2=Am2+z(n)*B{n};
end
vt=x0;

last_v1 = ones(max(size(A{1})),1) * inf;

while(norm(vt - last_v1, 2) > v_quadratic_error*norm(vt))
        last_v1 = vt;
             
        vt=(Am2)*D*(vt.*jump)+(sum((1-jump).*vt,1)*x0);
        
        
end

x=vt;


end





