function [bound_planes] = bound_planes_3d(x_n, msh_seg)
%
% function [bound_planes] = bound_planes_3d(x_n, msh_seg)
%
% Usage:   Computes the bounding planes to the triangular faces of the mesh elements. Mesh elements are tetrahedra
%
% Syntax:  [bound_planes] = bound_planes_3d(x_n, msh_seg)
%
% INPUT:
%    x_n        - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    msh_seg    - Segments of the given mesh | Required list of triangles [n_tri x 3]
%
% OUTPUT:
%    bound_planes   - planes bounding the msh_seg triangles containing
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

% Number of triangles
n_tri = size(msh_seg,1);

% Initialize carrier and side perpendicular planes' coefficients
% for all the triangles of the mesh | plane equation coefficients
tri_cplane   = zeros(size(msh_seg,1),4);
tri_p1_plane = zeros(size(msh_seg,1),4);
tri_p2_plane = zeros(size(msh_seg,1),4);
tri_p3_plane = zeros(size(msh_seg,1),4);

% Loop over the faces
for i = 1:n_tri
    
    % Face nodes: 1, 2, 3 -> A, B, C
    A = x_n(msh_seg(i,1),:); 
    B = x_n(msh_seg(i,2),:);
    C = x_n(msh_seg(i,3),:);
    
    % Face centroid
    centroid = (A+B+C) ./ 3;
    
    % Face edges
    ab = B - A;  ac = C - A;
    bc = C - B;  ca = A - C;
    
    %----------------------------------------------------------------------
    % Carrier plane : ax + by + cz + d = 0
    cp_norm = cross(ab,ac);
    cp_norm = cp_norm ./ norm(cp_norm);     % Make normal unit vector
    d_cp = - sum(cp_norm .* A);               
    
    % Store plane's coefficients
    tri_cplane(i,:) = [cp_norm d_cp];
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Perpendicular plane p1 normal vector
    p1_norm = cross(cp_norm,ab);
    p1_norm = p1_norm./norm(p1_norm);
    
    % Center of ab edge
    ab_cent = (B+A)./2;
    
    % Vector connecting ab center with face centroid
    abCC = centroid - ab_cent;
    
    % Check if dot product with p1 normal vector is positive
    if (sum(abCC.*p1_norm) < 0)
        % Flip p1 normal vector
        p1_norm = -p1_norm;
    end
    
    % Store p1 plane's coefficients
    d_p1 = -sum(p1_norm.*A);
    tri_p1_plane(i,:) = [p1_norm d_p1];
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    % Perpendicular plane p2 normal vector
    p2_norm = cross(cp_norm,bc);
    p2_norm = p2_norm./norm(p2_norm);
    
    
    % Center of bc edge
    bc_cent = (C+B)./2;
    
    % Vector connecting bc center with face centroid
    bcCC = centroid - bc_cent;
    
    % Check if dot product with p2 normal vector is positive
    if (sum(bcCC.*p2_norm) < 0)
        % Flip p2 normal
        p2_norm = -p2_norm;
    end

    % Store p2 plane's coefficients
    d_p2 = -sum(p2_norm.*B);
    tri_p2_plane(i,:) = [p2_norm d_p2];
    %----------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------- 
    % Perpendicular plane p3 normal vector 
    p3_norm = cross(cp_norm,ca);
    p3_norm = p3_norm./norm(p3_norm);
    
    % Center of ca edge
    ca_cent = (A+C)./2;
    
    % Vector connecting ca center with face centroid
    caCC = centroid - ca_cent;
    
    % Check if dot product with p3 normal vector is positive
    if (sum(caCC.*p3_norm) < 0)
        % Flip p3
        p3_norm = -p3_norm;
    end
    
    % Store p3 plane's coefficients
    d_p3 = -sum(p3_norm.*C);    
    tri_p3_plane(i,:) = [p3_norm d_p3];
    %----------------------------------------------------------------------  
 
end

% Construct bound_planes data structure for output
bound_planes.carrier = tri_cplane;
bound_planes.p1 = tri_p1_plane;
bound_planes.p2 = tri_p2_plane;
bound_planes.p3 = tri_p3_plane;


end

