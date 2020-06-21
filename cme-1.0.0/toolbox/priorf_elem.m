function [q, dq] = priorf_elem(x_n, x_s, opt_cme)
%
% function [q, dq] = priorf_elem(x_n, x_s, opt_cme)
%
% Usage:   Calculation of the CME prior weight function using R-functions.
%
% Syntax:  [q, dq] = priorf_elem(x_n, x_s, opt_cme)
%
% INPUT:
%    x_n      - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    x_s      - Sampling points | Format: [s x dims] where dims = 2 or 3
%    opt_cme  - Structure collecting the options for the CME approximants generation
%
% OUTPUT:
%    q        - Structure collecting the values of the prior weight functions for each sampling point x_s
%   dq        - Structure collecting the first gradient of the prior weight functions for each sampling point x_s
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

if (isfield(opt_cme,'enear'))
    enear   = opt_cme.enear;
else
    error('List of nodes enear affecting the cell is not defined')
end

if (isfield(opt_cme,'bd_segm_inner') && isfield(opt_cme,'bd_segm_all') ...
        && isfield(opt_cme,'bd_nodes') && isfield(opt_cme,'has_segm_outer'))
    bd_segm_inner  = opt_cme.bd_segm_inner;
    bd_segm_all    = opt_cme.bd_segm_all;
    has_segm_outer = opt_cme.has_segm_outer;
    bd_nodes       = opt_cme.bd_nodes;
    
else
    error('Auxiliary arrays at each node indicating the ring segments and boundary are not defined')
end

if (isfield(opt_cme,'spow') && isfield(opt_cme,'mpow'))
    spow     = opt_cme.spow;
    mpow     = opt_cme.mpow;
else
    error('spow and mpow are not defined')
end


if (isfield(opt_cme,'alpha') )
    alpha     = opt_cme.alpha;
else
    error('alpha parameter for 3D prior weights calculation are not defined')
end

if (isfield(opt_cme,'isConv'))
    isConv   = opt_cme.isConv;
else
    error('Flag indicating if the geometry is convex has not been defined.')
end


if (size(x_n,2) == 3)
    if (isfield(opt_cme,'tri_c') && isfield(opt_cme,'tri_p1') ...
            && isfield(opt_cme,'tri_p2') && isfield(opt_cme,'tri_p3'))
        
        tri_c   = opt_cme.tri_c;
        tri_p1  = opt_cme.tri_p1;
        tri_p2  = opt_cme.tri_p2;
        tri_p3  = opt_cme.tri_p3;
    else
        error('Carrier and perpendicular planes to 3D mesh triangle faces are required.')
    end
else
    tri_c = []; tri_p1 = []; tri_p2 = []; tri_p3 = [];
end


% Check partition of unity option
if (isfield(opt_cme,'pUnity'))
    pUnity   = opt_cme.pUnity;
else
    pUnity   = 1;
end

% Number of geometry nodes and dimension
[n_pts, n_dims] = size(x_n);

% Number of sample points in the element and dimension
[s_pts, s_dims] = size(x_s);     

if (n_pts < n_dims+1)
    error('The number of node points is bad set');
end

if (n_dims ~= s_dims)
    error('The node points and the samples have different dimension');
end

% Number of nearest neighbors affecting this element
enear_num = length(enear);


% We compute the prior
Zq  = zeros(s_pts,1);
q   = zeros(s_pts,enear_num);
dZq = zeros(s_pts,n_dims);
dq  = zeros(s_pts,enear_num,n_dims);
oneD= ones(1,n_dims);


% Loop over the nearest nodes affecting the element
for j=1:enear_num
    
    % Inner segments of the support domain
    segm_inner = bd_segm_inner{enear(j)};
  
    % If one of the ring edges of an interior node has a segment 
    % on the boundary OR domain is not convex
    if isConv==0 && has_segm_outer(j)==1 && bd_nodes(enear(j))==0
        
        % Compute prior function for all the segments too
        segm_all = bd_segm_all{enear(j)};
      
        % Prior function computed with Rfunction (equivalence)
        % spow-1 to preserve the "smoothness"
        [phi1,dphi1] = rfunction_equiv(segm_inner, x_n, x_s, spow-1, mpow, alpha, ...
                                       tri_c, tri_p1, tri_p2, tri_p3);
    
        [phi2,dphi2] = rfunction_equiv(segm_all, x_n, x_s, 1, mpow, alpha, ...
                                       tri_c, tri_p1, tri_p2, tri_p3);
      
        % Compute approximation function
        phi  = (phi1.*phi2);
    
        % Check for NaNs
        if (sum(isnan(phi)) > 0), disp('Nans in prior weight'), phi, pause; end
    
        % Compute first derivatives
        dphi = dphi1.*(phi2*oneD) + (phi1*oneD).*dphi2;
    
        % Check for NaNs
        if (sum(isnan(dphi(:))) > 0), disp('Nans in prior weight gradient'), dphi, dphi(isnan(dphi))=0; pause; end
      
    else  % For convex domains or for a full interior node OR for a boundary node in non convex
        
        [phi, dphi] = rfunction_equiv(segm_inner, x_n, x_s, spow, mpow, alpha, ...
                                      tri_c, tri_p1, tri_p2, tri_p3);
      
        if (sum(isinf(phi(:))) > 0), disp('Inf in prior weight'), phi, pause; end
    
        if (sum(isnan(dphi(:))) > 0), disp('Nans in prior weight gradient'), dphi, dphi(isnan(dphi))=0; pause; end
      
    end
    
    Zq        =  Zq +  phi;
    dZq       = dZq + dphi;
    q(:,j)    = phi;
    dq(:,j,:) = dphi;
end

% Prior functions are made a partition of unity using Shepard's method
if pUnity == 1
    aux = zeros(enear_num,n_dims);
    
    for k = 1:s_pts
        q(k,:)=q(k,:)/Zq(k);
        for kk = 1:n_dims
            aux(:,kk) = dq(k,:,kk)/Zq(k) - q(k,:)/Zq(k)*dZq(k,kk);
            if (sum(isnan(aux(:))) > 0), aux, pause; end
        end
        dq(k,:,:)=aux;
    end
end


end