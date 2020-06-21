function [] = check_cme_options(cme_options)
%
% function [] = check_cme_options(cme_options)
%
% Usage:   Check the CME options data structure for missing mandatory data
%
% Syntax:  [] = check_cme_options(cme_options)
%
% INPUT:
%    cme_options - Structure collecting the options for the cme approximants generation
%
% OUTPUT:
%    none
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


% Check the cme_options for errors
if (~isfield(cme_options, 'conn')) 
    error('Input connectivity list is required')
end

if (cme_options.dims == 2 && size(cme_options.conn,2) ~= 3)
    error('Input connectivity list not consistent with 2D geometry. Linear triangle cells required');
end

if (cme_options.dims == 3 && size(cme_options.conn,2) ~= 4)
    error('Input connectivity list not consistent with 3D geometry. Linear tetrahedral cells required');
end

if (~isfield(cme_options, 'nring') || cme_options.nring < 1 || cme_options.nring > 3)
    error('Number of support rings is required. Allowed Values [1 | 2 | 3]')
end


end

