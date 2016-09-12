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

function [x]=functionalPageRank_duplex(A,alpha,z)


% The following program takes in input a cell A of 2 sparse adjacency 
% matrices (ordered layers of the duplex networks A{1},A{2}) and computes 
%
%  - the Functional Multiplex PageRank x_i of node i given the vector of
%     influences z=[z^(1,0),z^(0,1),z^(1,1)]
%  
%   -alpha is the damping factor (usually taken 0.85)
%
% Options:
%
% The influences z=[z^(0,1),z^(1,0),z^(1,1)]
% Notes:
% The adjacency matrices A{1} and A{2} must have the same dimension and are
% unweighted
% For the  adjacency matrices we adopt the following convention
% A{1}(i,j)=1 if and only if j points to in layer 1 otherwise A{1}(i,j)=0
% A{2}(i,j)=1 if and only if j points to in layer 2 otherwise A{2}(i,j)=0
v_quadratic_error=0.0001;

B{1}=A{1}.*(1-A{2});
B{2}=(1-A{1}).*A{2};
B{3}=A{1}.*A{2};


x0=zeros(1,max(size(B{1}))); 

for jj=1:3
       
x0=x0+z(jj)*((sum(B{jj},1)>0)+(sum(B{jj},2)>0)');
end
x0=(x0>0);
norma=nnz(x0);
x0=x0/norma;
x0=x0';
D=zeros(1,max(size(B{1})));    
for jj=1:3,
D=D+sum(B{jj},1)*z(jj);
end
D=D+(D==0);
D=ones(1,max(size(B{1})))./D;
n=numel(D);
D=spdiags(D(:),0,n,n);


l=(sum(B{1},1)*z(1)+sum(B{2},1)*z(2)+sum(B{3},1)*z(3))>0;
jump=alpha*l';


vt=x0;

last_v1 = ones(max(size(A{1})),1) * inf;

while(norm(vt - last_v1, 2) > v_quadratic_error*norm(vt))
        last_v1 = vt;
             
        vt=(B{1}*z(1)+B{2}*z(2)+B{3}*z(3))*D*(vt.*jump)+(sum((1-jump).*vt,1)*x0);
        
        
end

x=vt;


end





