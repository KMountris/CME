%%
% Main program to test the Cell-based Maximum Entropy (CME) approximants
% in 2D and 3D geometries. CME can be used as an alternative to the 
% Moving Least Squares (MLS) and similar approximants for the construction 
% of weak-form meshless methods such as the Element Free Galerkin (EFG) method.
%
% The author of this program does not hold any responsibility for any
% inappropriate use and hopes this tool to be usefull in any way to the 
% broader scientific community. Please read the provided license before any use.
%
% author:  Konstantinos A. Mountris, PhD, University of Zaragoza, Zaragoza, Spain
% email:   konstantinos.mountris@gmail.com | kmountris@unizar.es
% web:     https://www.mountris.org
%
% 
% Last update: 20/10/2019
%
%
% References:
% [1] Mountris, KA, Bourantas, GC, Millán, D, et al. Cell‐based maximum entropy approximants for three‐dimensional domains: 
%     Application in large strain elastodynamics using the meshless total Lagrangian explicit dynamics method. 
%     Int J Numer Methods Eng. 2019; 1– 15. https://doi.org/10.1002/nme.6218 
%
% Disclaimer:
% This program has been based on the CME code for 2D problems that was previously implemented and distributed
% by Prof. Daniel Millán. Please seek information regarding the implementation of Prof. Daniel Millán and
% relative literature in the provided readme file.
%%

% Add paths to the workspace
resolve_cme_paths;

% Set the cell-based maximum entropy (cme) approximants options
opt_cme.verb  = 0;              % Verbose output, 0:off    1:on
opt_cme.inspect_support = 1;    % Inspect support domains, 0:off    1:on
opt_cme.dims  = 3;              % Dimensions, 2 or 3
opt_cme.grad  = 1;              % Computation of the Gradient, 0:OFF 1:ON
opt_cme.TolNR = 1.e-12;         % Newton-Raphson tolerance fo the nonlinear max-ent problem
opt_cme.isConv= 1;              % Flag indicating if the domain is convex:1 or non convex:0 (for priorfun computation)
opt_cme.nring = 2;              % Number ring of neighbors to consider for each node: 1, 2, 3.
opt_cme.spow  = 3;              % w^spow, where w is the approximation to the distance function (R-function)
opt_cme.mpow  = 2;              % Distance function derivative maximum degree of the approximation
opt_cme.alpha = 1.6;            % R-function zero set shape factor.

% Load mesh data
if (opt_cme.dims == 2)
    load('rectangle.mat', 'nodes', 'cells');
elseif (opt_cme.dims == 3)
    load('cube.mat', 'nodes', 'cells');
end

% Store the cells connectivity in the cme options.
opt_cme.conn = cells;

% Generate sampling points.
gp_per_cell = 1;                                                % For 2D case: 1, 3, 4, 7  |  For 3D case: 1, 4, 5
[gauss_p, gauss_w] = gauss_points(nodes, cells, gp_per_cell);

% Compute cme.
fprintf('Computing CME approximants...\n'); time = cputime;
cme = cme_approx(nodes, gauss_p, gp_per_cell, opt_cme);
fprintf('Elapsed time to compute CME approximants: %g s\n', cputime - time);

% Storing cme.
fprintf('Storing CME approximants in file...\n'); time = cputime;
store_cme_to_file(cme, 'cme_vals.txt', 'cme_neighs.txt');
fprintf('Elapsed time to store CME approximants: %g s\n', cputime - time);
