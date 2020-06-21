function [gam] = gamma_func_LS(x_n, x_s, near, q, lam0, dlam, t)
%
% function [gam] = gamma_func_LS(ndim, x_n, x_s, near, q, lam0, dlam, t)
%
% Usage:   Calculation of the max-ent basis with linear search.
%
% Syntax:  [gam] = gamma_func_LS(ndim, x_n, x_s, near, q, lam0, dlam, t)
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

dim = size(x_n,2);

lam = lam0 + t*dlam';
sum2=0;
for id=1:dim
    sum2=sum2 + lam(id)*(x_s(id)-x_n(near,id));
end

temp = q'.*exp(sum2);
Z = sum(temp);
gam = log(Z);

end