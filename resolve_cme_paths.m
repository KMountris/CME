function [] = resolve_cme_paths()
%
% function [] = resolve_cme_paths()
%
% Usage:   Adds the corresponding paths to the cme program in the Workspace
%
% Syntax:  [] = resolve_cme_paths()
%
% INPUT:
%    none
%
% OUTPUT:
%    none
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

str = which('resolve_cme_paths.m');
path = fileparts(str);
restoredefaultpath;

addpath(strcat(path,'/cme'));
addpath(strcat(path,'/cme/toolbox'));
addpath(strcat(path,'/data'));
addpath(strcat(path,'/gauss'));


end

