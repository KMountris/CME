function [D, dD] = rfunction_equiv(segments, x_n, s_p, spow, mpow, alpha, tri_c, tri_p1, tri_p2, tri_p3)
%
% function [D, dD] = rfunction_equiv(segments, x_n, s_p, spow, mpow, alpha, tri_c, tri_p1, tri_p2, tri_p3)
%
% Usage:    Evaluate the level set function using R-functions equivalence.
%
% Syntax:  [D, dD] = rfunction_equiv(segments, x_n, s_p, spow, mpow, alpha, tri_c, tri_p1, tri_p2, tri_p3)
%
% INPUT:
%    x_n    - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    x_s    - Sampling points | Format: [s x dims] where dims = 2 or 3
%    spow   - Smoothness modulation coefficient
%    mpow   - Normalization order coefficient
%    alpha  - R-function zero set shape factor
%    tri_c  - Carrier planes for tetrahedral cell faces
%    tri_p1 - 1st perpendicular plane to the carrier plane of a cell face
%    tri_p2 - 2nd perpendicular plane to the carrier plane of a cell face
%    tri_p3 - 3rd perpendicular plane to the carrier plane of a cell face
%
% OUTPUT:
%    D      - Level set function
%   dD      - Gradient of the level set function
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

% Number of sampling points and dimensions    
[sp_num, dims] = size(s_p);

% Initialize level set function and derivatives 
D    = zeros(sp_num,1);         
dD   = zeros(sp_num,dims);      

% Loop over number of sampling points
for j = 1:sp_num
  x_s = s_p(j,:);
  
  % Compute the non-negative edge weight and its derivative
  if (dims == 2)
      [rho, drho] = rhoweight_2D(segments, x_n, x_s);
  elseif (dims == 3)
      [rho, drho] = rhoweight_3D(segments, x_s, alpha, tri_c, tri_p1, tri_p2, tri_p3); 
  else 
      error('R function computation available only for 2D and 3D domains');
  end
  
  % Compute the R function and its gradient
  [w, dw]  = equivalence_formula(rho, drho, mpow);
  
  % Compute level set function and its gradient
  if (spow < 2)
      D(j)    = w;
      dD(j,:) = dw;
  else
      D(j)    = w^spow;
      dD(j,:) = spow*w^(spow-1)*dw;
  end
  
  % Check for NaN values at the vertices of the polygon
  if (sum(isnan(dD(j,:))) > 0) 
     dD(j,isnan(dD(j,:))) = 0;
  end

end

end


%%
%%
function [rho, drho] = rhoweight_2D(edges, x_n, x_s)
%
% function [rho, drho] = rhoweight_2D(edges, x_n, x_s)
%
% Usage:    Compute non-negative edge weight and its derivative for 2D polygon using triangle inequality
%
% Syntax:  [rho, drho] = rhoweight_2D(edges, x_n, x_s)
%
% INPUT:
%    edges  - Triangular cell edges
%    x_n    - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    x_s    - Sampling points | Format: [s x dims] where dims = 2 or 3
%
% OUTPUT:
%   rho     - Edge weight
%   drho    - Edge weight gradient
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

N    = size(edges,1);  % N edges

% Loop over the edges
for i = 1:N
  v1 = x_n(edges(i,1),:);     % 1st edge node 
  v2 = x_n(edges(i,2),:);     % 2nd edge node
  
  s1 = x_s - v1;
  s  = v2 - v1;
  d  = norm(s);     % length of edge
  id = 1/d;         % inverse length of edge
  
  xc = (v1+v2)*0.5; % edge center coords
  
  f  = (s1(1)*s(2) - s1(2)*s(1))*id;    % signed distance function
  df = [s(2), -s(1)]*id;                % and its gradient
  
  t  = ( (0.5*d)^2 - sum((x_s-xc).^2) )*id;   % trim function
  dt = -2*(x_s-xc)*id;                        % and its gradient
  

  val = sqrt(t*t+f^4);  
  rho(i)   = sqrt(f^2 + 0.25*(val-t)^2 );  % normalized function
  
  if abs(rho(i))<2*eps
    drho(i,:) = -df;
    rho(i) = 0;
  else
    drho(i,:) = (f*df + 0.25*(val-t) * ((t*dt+2*f^3*df)/val - dt) )/rho(i);
  end
  
end

end


%%
%%
function [rho, drho] = rhoweight_3D(faces, x_s, alpha, tri_c, tri_p1, tri_p2, tri_p3)
%
% function [rho, drho] = rhoweight_3D(faces, x_s, alpha, tri_c, tri_p1, tri_p2, tri_p3)
%
% Usage:    Compute non-negative face weight and its derivative for 
%           3D triangular faces extracted from tetrahedral cells
%
% Syntax:  [rho, drho] = rhoweight_3D(faces, x_s, alpha, tri_c, tri_p1, tri_p2, tri_p3)
%
% INPUT:
%    faces  - Tetrahedral cell faces
%    x_n    - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    x_s    - Sampling points | Format: [s x dims] where dims = 2 or 3
%    alpha  - R-function zero set shape factor
%    tri_c  - Carrier planes for tetrahedral cell faces
%    tri_p1 - 1st perpendicular plane to the carrier plane of a cell face
%    tri_p2 - 2nd perpendicular plane to the carrier plane of a cell face
%    tri_p3 - 3rd perpendicular plane to the carrier plane of a cell face
%
% OUTPUT:
%   rho     - Face weight
%   drho    - Face weight gradient
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

% Get the carrier plane for the given face
carrier_p = tri_c(faces(:,5),:);
f = carrier_p(:,1)*x_s(1) + carrier_p(:,2)*x_s(2) + carrier_p(:,3)*x_s(3) + carrier_p(:,4);
df = carrier_p(:,1:3);

% Get the perpendicular p1 plane for the given face
plane_p1 = tri_p1(faces(:,5),:);
p1 = plane_p1(:,1)*x_s(1) + plane_p1(:,2)*x_s(2) + plane_p1(:,3)*x_s(3) + plane_p1(:,4);
dp1 = plane_p1(:,1:3);

% Get the perpendicular p2 plane for the given face
plane_p2 = tri_p2(faces(:,5),:);
p2 = plane_p2(:,1)*x_s(1) + plane_p2(:,2)*x_s(2) + plane_p2(:,3)*x_s(3) + plane_p2(:,4);
dp2 = plane_p2(:,1:3);

% Get the perpendicular p3 plane for the given face
plane_p3 = tri_p3(faces(:,5),:);
p3 = plane_p3(:,1)*x_s(1) + plane_p3(:,2)*x_s(2) + plane_p3(:,3)*x_s(3) + plane_p3(:,4);
dp3 = plane_p3(:,1:3);

% Abbreviation to simplify conjuction notation
G = sqrt(p1.*p1 + p2.*p2 + alpha.*f.*f);
dG = (2*p1.*dp1 + 2*p2.*dp2 + alpha*2.*f.*df) ./ (2.*G);

% Trim volume
p12 = p1 + p2 - G;
p123_sq_root_part = sqrt(p12.*p12 + p3.*p3 + alpha.*f.*f);
t = p12 + p3 - p123_sq_root_part;

% Trim volume gradient
temp_val = 2*(p12).*(dp1+dp2-dG) + 2*p3.*dp3 + alpha*2.*f.*df;
dt = dp1+dp2+dp3-dG - (temp_val./(2 * p123_sq_root_part) );
    
% Normalized function rho
val1 = sqrt(t.*t + f.*f.*f.*f); 
normal_func = sqrt(f.*f + 0.25*(val1-t).*(val1-t) ); 
rho = normal_func'; 

% Gradient or normalized function
drho = (f.*df + 0.25*(val1-t) .* ((t.*dt + 2.*f.^3.*df)./val1 - dt) )./normal_func;

% Treatment for tiny values
tiny_val = abs(rho) < 2*eps;
if (nnz(tiny_val) ~= 0)
    drho(tiny_val) = -df(tiny_val);
end

end


%%
%%
function [w, dw] = equivalence_formula(rho, drho, mpow)
%
% function [w, dw] = equivalence_formula(rho, drho, mpow)
%
% Usage:    Compute R function equivalence and its gradient
%
% Syntax:  [w, dw] = equivalence_formula(rho, drho, mpow)
%
% INPUT:
%    faces  - weight
%    x_n    - weight gradient
%    mpow   - Normalization order coefficient
%
% OUTPUT:
%   w     - R function equivalence
%   dw    - R function equivalence gradient
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

ws  = sum(1./(rho.^mpow));
w   = ws^(-1/mpow);
dw  = ( rho.^(-mpow-1) * drho ) * w/ws;

end

