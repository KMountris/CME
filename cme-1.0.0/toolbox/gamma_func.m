function [gam, dgam, hgam, p_a, Z] = gamma_func(x_n, x_s, lam, near, q)
%
% function [gam, dgam, hgam, p_a, Z] = gamma_func(x_n, x_s, lam, near, q)
%
% Usage:   Calculation of the max-ent basis.
%
% Syntax:  [gam, dgam, hgam, p_a, Z] = gamma_func(x_n, x_s, lam, near, q)
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

sum2 = 0;
for id=1:dim
    sum2 = sum2 + lam(id)*(x_s(id)-x_n(near,id));
end

temp = q'.*exp(sum2);
Z = sum(temp);
p_a = temp/Z;
gam = log(Z);

dgam = zeros(1,dim);
for id=1:dim
    dgam(id) = sum( (x_s(id)-x_n(near,id)).*p_a(:) );
end

hgam = zeros(dim,dim);
for id = 1:dim
    for jd = 1:dim
        hgam(id,jd) = sum(p_a(:).*(x_s(id)-x_n(near,id)).*(x_s(jd)-x_n(near,jd)) ) - dgam(id)*dgam(jd);
    end
end

end
