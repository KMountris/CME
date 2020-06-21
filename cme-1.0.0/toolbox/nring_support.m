function [supports] = nring_support(x_n, opt_cme)
%
% function [supports] = nring_support(x_n, opt_cme)
%
% Usage:   Generates N-ring support domains for the geometry nodes x_n
%
% Syntax:  [supports] = nring_support(x_n, opt_cme)
%
% INPUT:
%    x_n        - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    opt_cme    - Structure collecting the options for the CME approximants generation
%
% OUTPUT:
%    supports   - Support domains data structure.
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
%%


% Plot the support domains interactively
if (isfield(opt_cme, 'inspect_support'))
    inspect_rings = opt_cme.inspect_support;          % 0 : no plot | 1 : plot supports 
else
    inspect_rings = 0;
end

% Get dimensions.
if (isfield(opt_cme, 'dims'))
    dims = opt_cme.dims;
else
    error('Domain dimensions have not been assigned')
end

% Get connectivity
if (isfield(opt_cme, 'conn'))
  conn = opt_cme.conn;
else
  error('Input cells connectivity has not been assigned')
end

% Get number of domain's rings.
if (isfield(opt_cme, 'nring'))
  nring = opt_cme.nring;
else
  nring = 2;
end

% Check if rings number is supported.
if (nring > 3)
  error('The maximum allowed number of rings is 3')
end

% Get domain's boundary segments and nodes indices
if (isfield(opt_cme, 'ebnd') && isfield(opt_cme, 'ibnd'))
    ebnd = opt_cme.ebnd;
    ibnd = opt_cme.ibnd;
else
    ebnd = freeBoundary(triangulation(conn,x_n));
    ibnd = unique(ebnd);
end

% Number of geometry nodes
n_num = size(x_n,1);

% Number of cells
el_num = size(conn,1);


% ------------- 1-RING ----------------------------------------------------
% Nodes occurence counter
n_count = zeros(n_num,1);

% Find node occurences in connectivity list
for el_id = 1:el_num
    % Increase count for nodes of current element
    n_count(conn(el_id,:)) = n_count(conn(el_id,:)) + 1;
end

% Initialize first ring elements for each node
first_ring = {[]};
for n_id = 1:n_num
    % first ring elements' ids zero for now
    first_ring{n_id} = zeros(1,n_count(n_id));
    n_count(n_id) = 1;
end

% Iterate over elements list
for el_id = 1:el_num
    % Get element's nodes
    el_nodes = conn(el_id,:);
    
    % Iterate over element's nodes
    for k = 1:length(el_nodes)
        
        % Current node's id
        n_id = el_nodes(k);
        
        % Position to store element id in the ring vector
        pos = n_count(n_id);
        
        % Store the element id
        first_ring{n_id}(pos) = el_id;
        
        % Increase the insert position pointer for the next element of the
        % ring
        n_count(n_id) = pos + 1;
    end
end
% -------------------------------------------------------------------------

% ------------- 2-RING and 3-RING------------------------------------------
% Find second elements' ring and third according 
% to what is required by the user
if (nring > 1)
    % Initialize second rings list
    second_ring = {[]};
    
    % Iterate over nodes
    for n_id = 1:n_num
        % Get neighbor nodes - all belonging to the first ring
        neigh_nodes = unique(conn(first_ring{n_id},:));
        % Second ring is the union of the neighbors first rings
        second_ring{n_id}= unique([first_ring{neigh_nodes}]);
    end
  
    if (nring == 2)
        last_ring = second_ring;
    else
        % Initialize third rings list
        third_ring = {[]};
        
        for n_id = 1:n_num
            % Get neighbor nodes - all belonging to the second ring now
            neigh_nodes = unique(conn(second_ring{n_id},:));
            % Third ring is the union of the second neighbors first rings
            third_ring{n_id} = unique([first_ring{neigh_nodes}]);
        end
        last_ring = third_ring;
    end
else
    % First and last rings are the same when nring is 1
    last_ring = first_ring;
end

%--------------------------------------------------------------------------

% Elements in last ring counter
el_count = zeros(el_num,1);
for j = 1:n_num
  el_count(last_ring{j}) = el_count(last_ring{j}) + 1;
end

% Map to each element, the nodes affecting it. That is the nodes for whom
% this element belongs in their n-ring support domain
nodes_in_elem = {[]};

% Iterate over elements
for el_id = 1:el_num
  nodes_in_elem{el_id} = zeros(1, el_count(el_id));
  el_count(el_id) = 1;
end

% Iterate over nodes
for j = 1:n_num
    % Elements in the jth node's last ring 
    jring = last_ring{j}; 
    
    % Iterate over the jring elements
    for k = 1:length(jring)
        el_id    = jring(k);
        jpos = el_count(el_id);
        nodes_in_elem{el_id}(jpos) = j;
        el_count(el_id) = jpos + 1;
    end
end

% List of boundary flags for the geometry nodes
bd_nodes = zeros(n_num,1);

% Cells for all the segments and only the interior segments 
%of each node's support domain
bd_segm_all   = {[]};
bd_segm_inner = {[]};
has_segm_outer = zeros(n_num,1);


% Initialize mesh boundary segments list: 3 edges - triangle / 4 faces - tetrahedron
% Format: segments number x (nodes of segment, belonging el. id)
msh_bnd = zeros( (dims+1)*el_num, dims+1 );

% Iterate over elements
for el_id = 1:el_num
    
    % Get all boundary segments (edges/faces) of the mesh
    if (dims == 2)
            % 1st edge n1-n2
            msh_bnd(3*(el_id-1)+1,1:dims) = sort([conn(el_id,1) conn(el_id,2)]);
            msh_bnd(3*(el_id-1)+1,dims+1) = el_id;
            
            % 2nd edge n2-n3
            msh_bnd(3*(el_id-1)+2,1:dims) = sort([conn(el_id,2) conn(el_id,3)]);
            msh_bnd(3*(el_id-1)+2,dims+1) = el_id;
            
            % 3rd edge n3-n1
            msh_bnd(3*(el_id-1)+3,1:dims) = sort([conn(el_id,3) conn(el_id,1)]);
            msh_bnd(3*(el_id-1)+3,dims+1) = el_id;
    else
            % 1st face n1-n2-n4
            msh_bnd(4*(el_id-1)+1,1:dims) = sort([conn(el_id,1) conn(el_id,2) conn(el_id,4)]);
            msh_bnd(4*(el_id-1)+1,dims+1) = el_id;
            
            % 2nd face n4-n2-n3
            msh_bnd(4*(el_id-1)+2,1:dims) = sort([conn(el_id,4) conn(el_id,2) conn(el_id,3)]);
            msh_bnd(4*(el_id-1)+2,dims+1) = el_id;
            
            % 3rd face n3-n2-n1
            msh_bnd(4*(el_id-1)+3,1:dims) = sort([conn(el_id,3) conn(el_id,2) conn(el_id,1)]);
            msh_bnd(4*(el_id-1)+3,dims+1) = el_id;
            
            % 4th face n1-n4-n3
            msh_bnd(4*(el_id-1)+4,1:dims) = sort([conn(el_id,1) conn(el_id,4) conn(el_id,3)]);
            msh_bnd(4*(el_id-1)+4,dims+1) = el_id;
    end
    
    
end

% Add boundary segments indices in the mesh boundary segments list
% Final format: segments number x (nodes of segment, belonging el. id, segment id) 
ids = 1:size(msh_bnd,1);
msh_bnd = [msh_bnd ids'];


% Compute bounding planes for all the boundary segments (triangles) in the
% 3d case
if (dims == 3)
    bound_planes = bound_planes_3d(x_n, msh_bnd);
end

% Sort nodes' indices in geometry boundary segments list for comparison with
% sorted msh_bnd later
ebnd = sort(ebnd,2);

% Iterate over nodes
for n_id = 1:n_num
    
    % Find all the boundary segments of the support domain
    [row_ids, ~] = find(msh_bnd(:,dims+1)==last_ring{n_id});
    sup_dom = msh_bnd(row_ids,:);
    
    % Find only the exterior boundary segments of the support domain
    [unV,IA,IB] = unique(sup_dom(:,1:dims),'rows');
    sup_dom_bnd = sup_dom(IA(histc(IB,1:size(unV,1)) == 1),:);
    bd_segm_all{n_id} = sup_dom_bnd;
    
    % Get only the inner segments (not on domain's boundary)
    [~,IA] = setdiff(sup_dom_bnd(:,1:dims),ebnd,'rows');
    support_bnd_inner = sup_dom_bnd(IA,:);
    bd_segm_inner{n_id} = support_bnd_inner;
   
    % Check if has boundary segments intersecting with the boundary of the domain
    if (size(bd_segm_all{n_id}, 1) ~= size(bd_segm_inner{n_id}, 1))
        has_segm_outer(n_id) = 1;
    end
    
    % Check if the node is on the geometry's boundary
    if any(n_id==ibnd), bd_nodes(n_id) = 1; end
    
    % View the support domain
    if (inspect_rings == 1)
        figure(2);clf
        if (dims == 2)
            hold on
            % Plot nodes and mesh
            plot(x_n(n_id,1),x_n(n_id,2),'m^','MarkerFaceColor','m','Markersize',12)
            trimesh(conn,x_n(:,1),x_n(:,2),zeros(n_num,1),'Facecolor','none','FaceLighting',...
                'none','EdgeColor','k','EdgeLighting','flat')
            
            % Plot first and last rings of support elements
            triplot(conn(last_ring{n_id}, :),x_n(:,1),x_n(:,2),'g-','LineWidth',2);
            triplot(conn(first_ring{n_id}, :),x_n(:,1),x_n(:,2),'b-','LineWidth',2);
            
            % Support domain boundary
            plot(x_n(support_bnd_inner(:,1:dims),1),x_n(support_bnd_inner(:,1:dims),2),'rs', ...
                'Markersize',7 ,'LineWidth',2);
            legend('eval. node','mesh','last ring','first ring','support border');
            title('Support domain inspection')
        else
            hold on    
            
            % Plot nodes
            view([-37.5 25]);
            plot3(x_n(n_id,1),x_n(n_id,2),x_n(n_id,3),'m^','MarkerFaceColor','m','Markersize',12)
            
            % Plot support domain segments
            patch('Vertices',x_n,'Faces',sup_dom_bnd(:,1:dims),...
                'FaceColor','g','FaceAlpha',0.8, 'linewidth',2);
            patch('Vertices',x_n,'Faces',support_bnd_inner(:,1:dims),...
                'FaceColor','r','FaceAlpha',0.8, 'linewidth',2);
            % Support domain inner segments' nodes
            plot3(x_n(support_bnd_inner(:,1:dims),1),x_n(support_bnd_inner(:,1:dims),2), ...
                x_n(support_bnd_inner(:,1:dims),3),'bs', 'Markersize',8,'MarkerFaceColor','b');
            % Plot mesh
            patch('Vertices',x_n,'Faces',get_tetra_faces(conn),...
                'FaceColor','none','FaceAlpha',0, 'linewidth',2);
            legend('eval.node','support seg','interior support seg','interior segm nodes','mesh');
            title('Support domain inspection')
            rotate3d on
        end
        axis equal
        pause
    end

end

% Construc supports data structure for output
supports.nodes_in_elem  = nodes_in_elem;
supports.bd_segm_all    = bd_segm_all;
supports.bd_segm_inner  = bd_segm_inner;
supports.has_segm_outer = has_segm_outer;
supports.bd_nodes       = bd_nodes;
supports.first_ring     = first_ring;
supports.last_ring      = last_ring;
if (dims == 3)
    supports.bound_planes = bound_planes;
end

end
