function [cme, opt_cme] = cme_approx(x_n, x_s, sp_in_cell, opt_cme)
%
% function [cme, opt_cme] = cme_approx(x_n, x_s, sp_per_el, opt_cme, type)
%
% Usage:   Generates Cell-based Maximum Entropy (CME) approximants for a set of sampling points
%
% Syntax:  [cme, opt_cme] = cme_approx(x_n, x_s, sp_per_el, opt_cme)
%
% INPUT:
%    x_n        - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    x_s        - Sampling points | Format: [s x dims] where dims = 2 or 3
%    sp_in_cell - Sampling points number in each cell | Format : integer value or [sp x 1]
%    opt_cme    - Structure collecting the options for the CME approximants generation
%
% OUTPUT:
%    cme        - The CME approximants for each sampling point (basis function; first derivative; nearest nodes to sampling points)
%    opt_cme    - Updated CME approximants generation options (containing additional data)
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

% Check and retrieve domain dimensions.
if (size(x_n,2) < 2 && size(x_n,2) > 3)
    error('Domain dimensions are not supported. Expected: [2 | 3]')
elseif (size(x_n,2) ~= size(x_s,2))
    error('Dimensions of geometry nodes (x_n) and sampling points (x_s) are not consistent')
else
    opt_cme.dims = size(x_n,2);
end

% Check for missing mandatory data from opt_cme.
check_cme_options(opt_cme);

% Find boundary segments and nodes indices if not available.
if (~isfield(opt_cme, 'ebnd') || ~isfield(opt_cme, 'ibnd'))
    ebnd = freeBoundary(triangulation(opt_cme.conn, x_n));
    
    % Boundary segments of the domain
    opt_cme.ebnd = ebnd;
    
    % Indices of domain boundary nodes
    opt_cme.ibnd = unique(ebnd);
end

% Construct n-ring support domains for geometry nodes.
supports = nring_support(x_n, opt_cme);

% Update cme options.
opt_cme.bd_segm_inner  = supports.bd_segm_inner;
opt_cme.bd_segm_all    = supports.bd_segm_all;
opt_cme.has_segm_outer = supports.has_segm_outer;
opt_cme.bd_nodes       = supports.bd_nodes;
opt_cme.first_ring     = supports.first_ring;
opt_cme.last_ring      = supports.last_ring;
opt_cme.nodes_in_elem  = supports.nodes_in_elem;

% Update specific cme options only for 3D case
if (size(x_n,2) == 3)
    opt_cme.tri_c     = supports.bound_planes.carrier;
    opt_cme.tri_p1    = supports.bound_planes.p1;
    opt_cme.tri_p2    = supports.bound_planes.p2;
    opt_cme.tri_p3    = supports.bound_planes.p3;
end

% Compute cme approximants.
cme = wrapper_cme(x_n, x_s, sp_in_cell, opt_cme);

end

