function [gp, gw] = gauss_points(x_n, cells, gp_per_cell)
%
% function [gp, gw] = gauss_points(x_n, cells, gp_per_cell)
%
% Usage:   Generates gauss points and their weights for 2D triangular or 3D tetrahedral cells
%
% Syntax:  [gp, gw] = gauss_points(x_n, cells, gp_per_cell)
%
% INPUT:
%    x_n         - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    cells       - Geometry cells  | Format: [c x dims] where dims = 3 or 4
%    gp_per_cell - Number of gauss points per each cell
%
% OUTPUT:
%    gp          - The coordinates of the gauss points
%    gw          - The weights of the gauss points
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

% Get nodes and cells dimensions
n_dim = size(x_n,2);
c_dim = size(cells,2);

if (n_dim == 2 && c_dim == 3)
    [gp, gw] = gauss_points_in_triangle(x_n, cells, gp_per_cell);
elseif (n_dim == 3 && c_dim == 4)
    [gp, gw] = gauss_points_in_tetrahedron(x_n, cells, gp_per_cell);
else
    error('Could not generate gauss points. Supported: 2D triangular cells or 3D tetrahedral cells');
end


end


%%
function [gp, gw] = gauss_points_in_triangle(x_n, cells, gp_per_cell)
%
% function [gp, gw] = gauss_points_in_triangle(x_n, cells, gp_per_cell)
%
% Usage:   Generates gauss points and their weights for 2D triangular cells
%
% Syntax:  [gp, gw] = gauss_points_in_triangle(x_n, cells, gp_per_cell)
%
% INPUT:
%    x_n         - Geometry nodes  | Format: [n x dims] where dims = 2
%    cells       - Geometry cells  | Format: [c x dims] where dims = 3
%    gp_per_cell - Number of gauss points per each cell
%
% OUTPUT:
%    gp          - The coordinates of the gauss points
%    gw          - The weights of the gauss points
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
%

% Check nodes dimension.
if size(x_n, 2) ~= 2, error('Nodes must be points in 2-dimensional space. Format: [nn x 2]'); end

% Check cells type.
if size(cells, 2) ~= 3, error('Connectivity of triangular cells expected. Format: [tt x 3]'); end

% Check requested gauss points per cell.
if (gp_per_cell ~= 1 && gp_per_cell ~= 3 && gp_per_cell ~= 4 && gp_per_cell ~= 7)
    error('Not supported number of gauss points per cell requested. Supported: 1, 3, 4, 7');
end


% Number of triangular cells.
trinum_tot = size(cells, 1);

% Generate the requested number of gauss points for each triangular cell
if (gp_per_cell == 1)
    g_num = trinum_tot;
    gp = zeros(g_num, 2); 
    gw = zeros(g_num, 1);
    
    for tt = 1:trinum_tot
        conn_tt = cells(tt,:);   
        x_tt = x_n(conn_tt,:);
        
        % Get the coordinates and weights of the gauss points
        gp(tt,:) = sum(x_tt)/3;
        gw(tt,:) = triangle_area(x_tt);
    end
    
elseif (gp_per_cell == 3)
    g_num = trinum_tot*3;
    gp = zeros(g_num,2); 
    gw = zeros(g_num,1); 
    
    for tt = 1:trinum_tot
        conn_tt = cells(tt,:);
        x_tt = x_n(conn_tt,:);
        
        % Get the coordinates of the gauss points
        gp((tt-1)*3+1,:) = (2/3)*x_tt(1,:) + (1/6)*x_tt(2,:) + (1/6)*x_tt(3,:);
        gp((tt-1)*3+2,:) = (1/6)*x_tt(1,:) + (2/3)*x_tt(2,:) + (1/6)*x_tt(3,:);
        gp((tt-1)*3+3,:) = (1/6)*x_tt(1,:) + (1/6)*x_tt(2,:) + (2/3)*x_tt(3,:);
        
        % Compute the area of the element.
        area_tt = triangle_area(x_tt);
        
        % Get the weights of the gauss points
        gw((tt-1)*3+1,:) = area_tt/3;
        gw((tt-1)*3+2,:) = area_tt/3;
        gw((tt-1)*3+3,:) = area_tt/3;
    end
    
elseif (gp_per_cell == 4)
        g_num = trinum_tot*4;
    gp = zeros(g_num,2);
    gw = zeros(g_num,1); 
    
    for tt = 1:trinum_tot     
        conn_tt = cells(tt,:);
        x_tt = x_n(conn_tt,:);
        
        % Get the coordinates of the gauss points
        gp((tt-1)*4+1,:) = 0.6*x_tt(1,:) + 0.2*x_tt(2,:) + 0.2*x_tt(3,:);
        gp((tt-1)*4+2,:) = 0.2*x_tt(1,:) + 0.6*x_tt(2,:) + 0.2*x_tt(3,:);
        gp((tt-1)*4+3,:) = 0.2*x_tt(1,:) + 0.2*x_tt(2,:) + 0.6*x_tt(3,:);
        gp((tt-1)*4+4,:) = (1/3)*x_tt(1,:) + (1/3)*x_tt(2,:) + (1/3)*x_tt(3,:);
        
        % Compute the area of the element.
        area_tt = triangle_area(x_tt);
        
        % Get the weights of the gauss points
        gw((tt-1)*4+1,:) = 0.52083333333333*area_tt;
        gw((tt-1)*4+2,:) = 0.52083333333333*area_tt;
        gw((tt-1)*4+3,:) = 0.52083333333333*area_tt;
        gw((tt-1)*4+4,:) = -0.5625*area_tt;
    end  
    
elseif (gp_per_cell == 7)
    g_num = trinum_tot*7;
    gp = zeros(g_num,2);
    gw = zeros(g_num,1); 
    
    for tt = 1:trinum_tot     
        conn_tt = cells(tt,:);
        x_tt = x_n(conn_tt,:);  
        
        % Get the coordinates of the gauss points
        gp((tt-1)*7+1,:) = x_tt(3,:);
        gp((tt-1)*7+2,:) = 0.5*x_tt(1,:) + 0.5*x_tt(3,:);
        gp((tt-1)*7+3,:) = x_tt(1,:);
        gp((tt-1)*7+4,:) = 0.5*x_tt(1,:) + 0.5*x_tt(2,:);
        gp((tt-1)*7+5,:) = x_tt(2,:);
        gp((tt-1)*7+6,:) = 0.5*x_tt(2,:) + 0.5*x_tt(3,:);
        gp((tt-1)*7+7,:) = (1/3)*x_tt(1,:) + (1/3)*x_tt(2,:) + (1/3)*x_tt(3,:);
        
        % Compute the area of the element.
        area_tt = triangle_area(x_tt);
        
        % Get the weights of the gauss points
        gw((tt-1)*7+1,:) = 0.05*area_tt;
        gw((tt-1)*7+2,:) = 0.133333333333333*area_tt;
        gw((tt-1)*7+3,:) = 0.05*area_tt;
        gw((tt-1)*7+4,:) = 0.133333333333333*area_tt;
        gw((tt-1)*7+5,:) = 0.05*area_tt;
        gw((tt-1)*7+6,:) = 0.133333333333333*area_tt;
        gw((tt-1)*7+7,:) = 0.45*area_tt;
    end
    
end


end


%%
function [gp, gw] = gauss_points_in_tetrahedron(x_n, cells, gp_per_cell)
%
% function [gp, gw] = gauss_points_in_tetrahedron(x_n, cells, gp_per_cell)
%
% Usage:   Generates gauss points and their weights for 3D tetrahedral cells
%
% Syntax:  [gp, gw] = gauss_points_in_tetrahedron(x_n, cells, gp_per_cell)
%
% INPUT:
%    x_n         - Geometry nodes  | Format: [n x dims] where dims = 3
%    cells       - Geometry cells  | Format: [c x dims] where dims = 4
%    gp_per_cell - Number of gauss points per each cell
%
% OUTPUT:
%    gp          - The coordinates of the gauss points
%    gw          - The weights of the gauss points
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
%

% Check nodes dimension.
if size(x_n, 2) ~= 3, error('Nodes must be points in 3-dimensional space. Format: [nn x 3]'); end

% Check elements type.
if size(cells, 2) ~= 4, error('Connectivity of tetrahedral elements expected. Format: [tt x 4]'); end

% Check requested gauss points per cell.
if (gp_per_cell ~= 1 && gp_per_cell ~= 4 && gp_per_cell ~= 5)
    error('Not supported number of gauss points per cell requested. Supported: 1, 4, 5');
end


% Number of tetrahedral cells.
tetnum_tot = size(cells, 1);

if (gp_per_cell == 1)
    g_num = tetnum_tot;
    gp = zeros(g_num, 3);    
    gw = zeros(g_num, 1);  
    
    for tt = 1:tetnum_tot
        conn_tt = cells(tt,:);
        x_tt = x_n(conn_tt,:);
            
        % Get the coordinates and weights of the gauss points
        gp(tt,:) = 1/4*x_tt(1,:) + 1/4*x_tt(2,:) + 1/4*x_tt(3,:) + 1/4*x_tt(4,:);
        gw(tt,:) = tetrahedron_volume(x_tt);
    end

elseif (gp_per_cell == 4) 
    g_num = tetnum_tot*4;
    gp = zeros(g_num, 3); 
    gw = zeros(g_num, 1); 

    for tt = 1:tetnum_tot
        conn_tt = cells(tt,:);
        x_tt = x_n(conn_tt,:);
           
        % Get the coordinates of the gauss points
        alpha = 0.58541020; 
        beta = 0.13819660;
        gp((tt-1)*4 + 1,:) = alpha*x_tt(1,:) + beta*x_tt(2,:) + beta*x_tt(3,:) + beta*x_tt(4,:);
        gp((tt-1)*4 + 2,:) = beta*x_tt(1,:) + alpha*x_tt(2,:) + beta*x_tt(3,:) + beta*x_tt(4,:);
        gp((tt-1)*4 + 3,:) = beta*x_tt(1,:) + beta*x_tt(2,:) + alpha*x_tt(3,:) + beta*x_tt(4,:);
        gp((tt-1)*4 + 4,:) = beta*x_tt(1,:) + beta*x_tt(2,:) + beta*x_tt(3,:) + alpha*x_tt(4,:);
        
        % Compute tetrahedron volume
        volume_el = tetrahedron_volume(x_tt);
        
        % Get the weights of the gauss points
        gw((tt-1)*4 + 1,:) = volume_el/4;
        gw((tt-1)*4 + 2,:) = volume_el/4;
        gw((tt-1)*4 + 3,:) = volume_el/4;
        gw((tt-1)*4 + 4,:) = volume_el/4;

    end

elseif (gp_per_cell == 5) 
    g_num = tetnum_tot*5;
    gp = zeros(g_num, 3); 
    gw = zeros(g_num, 1); 
 
    for tt = 1:tetnum_tot
        conn_tt = cells(tt,:);
        x_tt = x_n(conn_tt,:);
           
        % Get the coordinates of the gauss points
        gp((tt-1)*5 + 1,:) = 1/4*x_tt(1,:) + 1/4*x_tt(2,:) + 1/4*x_tt(3,:) + 1/4*x_tt(4,:);
        gp((tt-1)*5 + 2,:) = 1/2*x_tt(1,:) + 1/6*x_tt(2,:) + 1/6*x_tt(3,:) + 1/6*x_tt(4,:);
        gp((tt-1)*5 + 3,:) = 1/6*x_tt(1,:) + 1/2*x_tt(2,:) + 1/6*x_tt(3,:) + 1/6*x_tt(4,:);
        gp((tt-1)*5 + 4,:) = 1/6*x_tt(1,:) + 1/6*x_tt(2,:) + 1/2*x_tt(3,:) + 1/6*x_tt(4,:);
        gp((tt-1)*5 + 5,:) = 1/6*x_tt(1,:) + 1/6*x_tt(2,:) + 1/6*x_tt(3,:) + 1/2*x_tt(4,:);
        
        % Compute tetrahedron volume
        volume_el = tetrahedron_volume(x_tt);
        
        % Get the weights of the gauss points
        gw((tt-1)*5 + 1,:) = volume_el*(-4/5);
        gw((tt-1)*5 + 2,:) = volume_el*(9/20);
        gw((tt-1)*5 + 3,:) = volume_el*(9/20);
        gw((tt-1)*5 + 4,:) = volume_el*(9/20);
        gw((tt-1)*5 + 5,:) = volume_el*(9/20);
    end       
end

end


%%
function [a] = triangle_area(x)
%
% function [a] = triangle_area(x)
%
% Usage:   Computes the area of a 2D triangular cell
%
% Syntax:  [a] = triangle_area(x)
%
% INPUT:
%    x   - Coordinates of the triangular cell's nodes  | Format: [3 x 2]
%
% OUTPUT:
%    a   - The triangular cell area
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
%

% Check triangle nodes coordinates dimensions.
if (size(x,1) ~= 3 && size(x, 2) ~= 2)
    error('Cannot compute area of given triangular cell. Expected format: [3, 2]'); 
end

% Computation of triangular cell area.
a = x(2,1)*x(3,2)-x(3,1)*x(2,2) + ...
    x(3,1)*x(1,2)-x(1,1)*x(3,2) + ...
    x(1,1)*x(2,2)-x(2,1)*x(1,2);
 
a = 0.5*abs(a);

end


%%
function [v] = tetrahedron_volume(x)
%
% function [v] = tetrahedron_volume(x)
%
% Usage:   Computes the volume of a 3D tetrahedral cell
%
% Syntax:  [v] = tetrahedron_volume(x)
%
% INPUT:
%    x   - Coordinates of the tetrahedral cell's nodes  | Format: [4 x 3]
%
% OUTPUT:
%    v   - The tetrahedral cell volume
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
%

% Check tetrahedron nodes coordinates dimensions.
if (size(x,1) ~= 4 && size(x, 2) ~= 3)
    error('Cannot compute volume of given tetrahedral cell. Expected format: [4, 3]'); 
end

x1 = x(1,1); y1 = x(1,2); z1 = x(1,3);
x2 = x(2,1); y2 = x(2,2); z2 = x(2,3);
x3 = x(3,1); y3 = x(3,2); z3 = x(3,3);
x4 = x(4,1); y4 = x(4,2); z4 = x(4,3);

a1 = y1*(z2-z3)-z1*(y2-y3)+(y2*z3-y3*z2);
a2 = x1*(z2-z3)-z1*(x2-x3)+(z3*x2-z2*x3);
a3 = x1*(y2-y3)-y1*(x2-x3)+(x2*y3-x3*y2);
a4 = x1*(y2*z3-y3*z2)-y1*(x2*z3-x3*z2)+z1*(x2*y3-x3*y2);

v = abs((y4*a2-z4*a3+a4-x4*a1)/6);


end



