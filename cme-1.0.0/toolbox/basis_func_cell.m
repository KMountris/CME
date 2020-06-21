function [p, dp] = basis_func_cell(x_n, x_s, opt_cme)
%
% function [p, dp] = basis_func_cell(x_n, x_s, opt_cme)
%
% Usage:   Calculation of the CME basis function using the Newton-Rhapson method.
%
% Syntax:  [p, dp] = basis_func_cell(x_n, x_s, opt_cme)
%
% INPUT:
%    x_n      - Geometry nodes  | Format: [n x dims] where dims = 2 or 3
%    x_s      - Sampling points | Format: [s x dims] where dims = 2 or 3
%    opt_cme  - Structure collecting the options for the CME approximants generation
%
% OUTPUT:
%    p        - Structure collecting the values of the basis functions for each sampling point x_s
%   dp        - Structure collecting the first gradient of the basis functions for each sampling point x_s
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

warning('off')

if isfield(opt_cme,'TolNR')
  TolNR   = opt_cme.TolNR;
else
  error('Tolerance TolNR for Newton-Raphson algorithm is not defined')
end

if isfield(opt_cme,'enear')
  enear   = opt_cme.enear;
else
  error('List of nodes enear affecting the cell is not defined')
end

[nPts,nDim] = size(x_n);
[sPts,sDim] = size(x_s);     %number of sampling points in the cell
if nPts<nDim+1
  error('The number of node points is set wrong');
end
if nDim ~= sDim
  error('The node points and the sampling points have different dimension');
end

eNN      = length(enear);   %number of nearest neighbors affecting this cell
i_fail   = zeros(sPts,1);
max_iter = 100;

% Basis functions initialization.
p   = zeros(sPts,eNN);
dp  = zeros(sPts,eNN,nDim);
aux = zeros(eNN,nDim);

% Prior computation using Rfunctions 
[q,dq] = priorf_elem(x_n,x_s,opt_cme);

% Cell-based max-ent functions computation
for k = 1:sPts
  
  % Newton-Rhapson (N-R) method
  x  = x_s(k,:);
  lam= zeros(1,nDim);
  R  = 10*ones(1,nDim);
  niter=0;
  
  % N-R Regular iteration
  while (norm(R)>TolNR) && niter <= max_iter
    [~,R,J,p_a]=gamma_func(x_n,x,lam,enear,q(k,:));
    if (abs(rcond(J)))<1e-8
      p(k,:)=NaN([eNN 1]);
      i_fail(k)=1;
      break
    end
    dlam=-J\R';
    if sum(dlam.*R')>0
      p(k,:)=NaN([eNN 1]);
      i_fail(k)=1;
      break
    end
    lam=lam+dlam';
    niter=niter+1;
  end
  if niter >=max_iter
    p(k,:)=NaN([eNN 1]);
    i_fail(k)=1;
  end
  if any(isnan(p_a))
    p(k,:)=NaN([eNN 1]);
    i_fail(k)=1;
  end
  
  % N-R with Line-Search when things go wrong
  if i_fail(k) == 1
    i_fail(k)=0;
    lam=zeros(1,nDim);
    R=10*ones(1,nDim);
    niter=0;
    while (norm(R)>TolNR) && niter <= max_iter
      [~,R,J,p_a]=gamma_func(x_n,x,lam,enear,q(k,:));
      dlam=-J\R';
      
      %Line-Search
      if sum(dlam.*R')>0
        dlam=-dlam;
      end
      opts=optimset('TolX',1e-8,'MaxIter',1000);
      t = fminbnd(@(t) gamma_func_LS(x_n,x,enear,q(k,:),lam,dlam,t),0,1,opts);
      lam=lam+t*dlam';
      niter=niter+1;
    end
    
    if niter >=max_iter
      fprintf(1,'Newton-Rhapson Line-Search Failed, max num iterations=%d  |R|=%12.2e  TolNR=%12.2e\n',niter,norm(R),TolNR);
      p(k,:)=NaN([eNN 1]);
      i_fail(k)=1;
      
      figure(50);clf
      hold on
      if nDim==1
        plot(x_n,0*x_n,'ro')
        plot(x_n(enear),0*x_n(enear),'g*',x,0,'mx','Markersize',10 ,'LineWidth',2);
      elseif nDim==2
        plot(x_n(:,1),x_n(:,2),'ro')
        plot(x_n(enear,1),x_n(enear,2),'g*',x(1),x(2),'mx','Markersize',10 ,'LineWidth',2);
        trimesh(opt_cme.conn,x_n(:,1),x_n(:,2),zeros(size(x_n,1),1),'Facecolor','none','FaceLighting',...
                'flat','EdgeColor','k','EdgeLighting','flat');
      else
        plot3(x_n(:,1),x_n(:,2),x_n(:,3),'ro')
        plot3(x_n(enear,1),x_n(enear,2),x_n(enear,3),'g*',x(1),x(2),x(3),'mx','Markersize',10 ,'LineWidth',2);
      end
      hold off
      axis equal
      pause(0.1)
    end
    
    if any(isnan(p_a))
      disp('Newton-Rhapson Line-Search Failed, NaNs')
      p(k,:)=NaN([eNN 1]);
      p(k,:)
      i_fail(k)=1;

      figure(60);clf
      hold on
      if nDim==1
        plot(x_n,0*x_n,'ro')
        plot(x_n(enear),0*x_n(enear),'g*',x,0*x,'cx','Markersize',10 ,'LineWidth',2);
      elseif nDim==2
        plot(x_n(:,1),x_n(:,2),'ro')
        plot(x_n(enear,1),x_n(enear,2),'g*',x(1),x(2),'cx','Markersize',10 ,'LineWidth',2);
      else
        plot3(x_n(:,1),x_n(:,2),x_n(:,3),'ro')
        plot3(x_n(enear,1),x_n(enear,2),x_n(enear,3),'g*',x(1),x(2),x(3),'mx','Markersize',10 ,'LineWidth',2);
      end
      hold off
      axis equal
      pause(0.1)
    end
    
  end
  
  %Checking for fails
  if i_fail(k) == 0
    p(k,:)=p_a;
    x_m_xa=zeros(nDim,eNN);
    for kk = 1:nDim
      x_m_xa(kk,:) = x(kk)-x_n(enear,kk);
    end
    aux0 = exp(lam*x_m_xa)';
    aux0 = aux0/sum(q(k,:)'.*aux0);
    dr_dx = zeros(nDim,nDim);
    for kk = 1:nDim
      for jj = 1:nDim
        dr_dx(kk,jj)= sum(aux0.*dq(k,:,jj)'.*x_m_xa(kk,:)');
      end
    end
    dr_dx=dr_dx+eye(nDim);
    dlamdx = -inv(J)*dr_dx ;
    aux1 = dlamdx'*x_m_xa;
    for kk = 1:nDim
      aux2 = sum(aux0.*dq(k,:,kk)');
      aux(:,kk) = aux0.*dq(k,:,kk)' + p_a.*aux1(kk,:)' - p_a * aux2;
    end
    dp(k,:,:) = aux;
  end
  
  
end

