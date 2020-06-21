function [faces] = get_tetra_faces(cells)
%
% function [faces] = get_tetra_faces(cells)
%
% Usage:   Extract the triangular faces of a tetrahedral cell.
%
% Syntax:  [faces] = get_tetra_faces(cells)
%
%
% Author: Konstantinos A. Mountris, PhD, University of Zaragoza, Zaragoza, Spain
% email:  konstantinos.mountris@gmail.com | kmountris@unizar.es
% web:    https://www.mountris.org
% 
% Last update: 20/10/2019
%
%
% References:
% [1] Mountris, KA, Bourantas, GC, Millán, D, et al. Cell‐based maximum entropy approximants for three‐dimensional domains: 
%     Application in large strain elastodynamics using the meshless total Lagrangian explicit dynamics method. 
%     Int J Numer Methods Eng. 2019; 1– 15. https://doi.org/10.1002/nme.6218 
%
%%

fc = [1 2 4
      1 3 2
      2 3 4
      1 4 3];
    
cells_num = size(cells,1);    
faces = zeros(4*cells_num, 3);    
pos = 1;

for kk = 1:cells_num
    cl = cells(kk,:);
    faces(pos:pos+3,:) = cl(fc);
    pos = pos + 4;
end

end
