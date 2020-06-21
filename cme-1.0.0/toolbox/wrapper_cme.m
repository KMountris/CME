function [cme] = wrapper_cme(x_n, x_s, sp_in_cell, opt_cme)
%
% function [cme] = wrapper_cme(x_n, x_s, sp_per_el, opt_cme)
%
% Usage:   Computes the CME approximants for the sampling points (x_s)
%
% Syntax:  [cme] = wrapper_cme(x_n, x_s, sp_per_el, opt_cme)
%
% INPUT:
%    x_n        - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    x_s        - Sampling points | Format: [s x dims] where dims = 2 or 3
%    sp_in_cell - Sampling points number in each cell | Format : integer value or [sp x 1]
%    opt_cme    - Structure collecting the options for the CME approximants generation
%
% OUTPUT:
%    cme        - The CME approximants for each sampling point (basis function; first derivative; nearest nodes to sampling points)
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


if (nargout > 1)
    error('The number of output is too big');
end

if (isfield(opt_cme,'nodes_in_elem'))
    nodes_in_elem = opt_cme.nodes_in_elem;
else
    error('Auxiliary arrays at each node indicating "nodes_in_elem" is not defined')
end

if (isfield(opt_cme, 'conn'))
    conn = opt_cme.conn;
else
    error('Input connectivity list is not defined')
end

% Number of sampling points
s_pts   = size(x_s,1);

% Geometry dimensions
dims   = opt_cme.dims;      

% Number of connectivity's elements
elem_num = size(conn,1);

% Elements indices
elem_ids = (1:elem_num)';

% Sample points indices
sp_ids = (1:size(x_s,1))';

% Index of the last sample point per element
if (length(sp_in_cell) == 1)
    last_sp_id_in_el = sp_in_cell*elem_ids;
elseif (length(sp_in_cell) == s_pts)
    last_sp_id_in_el = cumsum(sp_in_cell);
else
    error('x_s and sp_in_cell have different size')
end

% === Compute basis functions in each element containing sample points ===

% Initialize approximation function and 1st derivative values on the
% sample points
p_cme  = {zeros(s_pts,1)};
dp_cme = {zeros(s_pts,1)};

% Initialize nearest neighbors nodes to sample points
s_near = {zeros(s_pts,1)};

% Initialize id of the first sample point id in an element
first_sp_id = 0;

% Iterate over elements
for e = 1 : elem_num
     % Global element index
     iel = elem_ids(e);
     
     % Nearest nodes to the e-th element
     enear = nodes_in_elem{iel};
     opt_cme.enear = enear;    
  
     % Number of sample points in the e-th triangle
     kPts = last_sp_id_in_el(e) - first_sp_id;   
     
     % Indices of the sample points belonging to the e-th triangle
     id_samp = (first_sp_id+1):(first_sp_id+kPts);
     id_samp = sp_ids(id_samp);
     
     % Sample points coordinates
     xe_s = x_s(id_samp,:);   
     
     % cell-based max-ent functions and 1st derivatives
     [p, dp] = basis_func_cell(x_n, xe_s, opt_cme);
 
     if (opt_cme.verb == 1)
         fprintf('computed cme for cell: %d\n', e);
     end
     
     for j=1:kPts
         s = id_samp(j);  
         
         p_cme{s} = p(j,:)';
         if (dims == 2)
             dp_cme{s} = ([dp(j,:,1);dp(j,:,2)]');
         else
             dp_cme{s} = ([dp(j,:,1);dp(j,:,2);dp(j,:,3)]');
         end
         s_near{s} = enear;
     end
     first_sp_id = last_sp_id_in_el(e);
end

cme.phi  = p_cme;
cme.dphi = dp_cme;
cme.neighs  = s_near;
