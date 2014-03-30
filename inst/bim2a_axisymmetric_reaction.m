## Copyright (C) 2006-2014  Carlo de Falco, Massimiliano Culpo
##
## This file is part of:
##     BIM - Diffusion Advection Reaction PDE Solver
##
##  BIM is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  BIM is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with BIM; If not, see <http://www.gnu.org/licenses/>.
##
##  author: Carlo de Falco     <cdf _AT_ users.sourceforge.net>
##  author: Massimiliano Culpo <culpo _AT_ users.sourceforge.net>
##  author: Matteo porro       <meoo85 _AT_ users.sourceforge.net>
##  author: Emanuela Abbate    <emanuela.abbate _AT_ mail.polimi.it>

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{C}]} = @
## bim2a_axisymmetric_reaction(@var{mesh},@var{delta},@var{zeta}) 
##
## Build the lumped finite element mass matrix for a diffusion
## problem in cylindrical coordinates with axisymmetric configuration.
##
## The equation taken into account is:
##
## @var{delta} * @var{zeta} * u = f
## 
## where @var{delta} is an element-wise constant scalar function, while
## @var{zeta} is a piecewise linear conforming scalar function.
##
## @seealso{bim2a_rhs, bim2a_axisymmetric_advection_diffusion, 
## bim2a_axisymmetric_laplacian, bim2a_reaction, bim1a_reaction, bim3a_reaction}
## @end deftypefn

function [C] = bim2a_axisymmetric_reaction(mesh,delta,zeta)

  ## Check input
  if nargin != 3
    error("bim2a_axisymmetric_reaction: wrong number of input parameters.");
  elseif !(isstruct(mesh)     && isfield(mesh,"p") &&
	   isfield (mesh,"t") && isfield(mesh,"e"))
    error("bim2a_axisymmetric_reaction: first input is not a valid mesh structure.");
  elseif !(all(mesh.p(1,:) >= 0) || all(mesh.p(1,:) <= 0))
    error("bim2a_axisymmetric_reaction: the input mesh cannot intersect the rotation axis r=0.");
  endif
  
  nnodes = size(mesh.p,2);
  nelem  = size(mesh.t,2);

  r = abs (mesh.p(1,:));

  ## Turn scalar input to a vector of appropriate size
  if isscalar(delta)
    delta = delta*ones(nelem,1);
  endif
  if isscalar(zeta)
    zeta = zeta*ones(nnodes,1);
  endif

  if !( isvector(delta) && isvector(zeta) )
    error("bim2a_axisymmetric_reaction: coefficients are not valid vectors.");
  elseif length(delta) != nelem
    error("bim2a_axisymmetric_reaction: length of alpha is not equal to the number of elements.");
  elseif length(zeta) != nnodes
    error("bim2a_axisymmetric_reaction: length of gamma is not equal to the number of nodes.");
  endif

  wjacdet   = mesh.wjacdet(:,:);
  coeff     = zeta(mesh.t(1:3,:));
  coeffe    = delta(:);
  
  ## Local matrix	
  Blocmat = zeros(3,nelem);	
  for inode = 1:3
    Blocmat(inode,:) = coeffe'.*coeff(inode,:).*wjacdet(inode,:) .* r(mesh.t(inode,:));
  endfor
  
  gnode = (mesh.t(1:3,:));
  
  ## Global matrix
  C = sparse(gnode(:),gnode(:),Blocmat(:));

endfunction

%!shared mesh,delta,zeta,nnodes,nelem
% x = y = linspace(0,1,4);
% [mesh] = msh2m_structured_mesh(x,y,1,1:4);
% [mesh] = bim2c_mesh_properties(mesh);
% nnodes = columns(mesh.p);
% nelem  = columns(mesh.t);
% delta  = ones(columns(mesh.t),1);
% zeta   = ones(columns(mesh.p),1);
%!test
% [C] = bim2a_axisymmetric_reaction(mesh,delta,zeta);
% assert(size(C),[nnodes, nnodes]);
%!test
% [C1] = bim2a_axisymmetric_reaction(mesh,3*delta,zeta);
% [C2] = bim2a_axisymmetric_reaction(mesh,delta,3*zeta);
% assert(C1,C2);
%!test
% [C1] = bim2a_axisymmetric_reaction(mesh,3*delta,zeta);
% [C2] = bim2a_axisymmetric_reaction(mesh,3,1);
% assert(C1,C2);

%!test
%! n = 20;
%! [mesh] = msh2m_structured_mesh(linspace(1,2,n+1),linspace(0,1,n+1),1,1:4);
%! mesh   = bim2c_mesh_properties(mesh);
%! uex    = @(r,z) exp(r) .* exp(1-z);
%! duexdr = @(r,z) uex(r,z);
%! d2uexdr2 = @(r,z) uex(r,z);
%! duexdz = @(r,z) -uex(r,z);
%! d2uexdz2 = @(r,z) uex(r,z);
%! Dnodes = bim2c_unknowns_on_side(mesh,[1,2,3,4]);
%! Nnodes    = columns(mesh.p);
%! Nelements = columns(mesh.t);
%! Varnodes  = setdiff(1:Nnodes,Dnodes);
%! D = 1; vr = 1; vz = 1; sigma = 1;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! delta  = sigma*ones(columns(mesh.t),1);
%! zeta   = ones(columns(mesh.p),1);
%! beta   = 1/D*[vr*ones(1,Nelements); vz*ones(1,Nelements)];
%! f = @(r,z) -D./r.*duexdr(r,z) - D.*d2uexdr2(r,z) ...
%!           + vr./r .* uex(r,z) + vr * duexdr(r,z) ...
%!           - D.*d2uexdz2(r,z) + vz * duexdz(r,z) ...
%!           + sigma * uex(r,z);
%! rhs    = bim2a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh.p(1,:), mesh.p(2,:)));
%! S = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! C = bim2a_axisymmetric_reaction(mesh,delta,zeta);
%! S += C;
%! u = zeros(Nnodes,1); u(Dnodes) = uex(mesh.p(1,Dnodes), mesh.p(2,Dnodes));
%! u(Varnodes) = S(Varnodes,Varnodes)\(rhs(Varnodes) - S(Varnodes,Dnodes)*u(Dnodes));
%! assert(u,uex(mesh.p(1,:), mesh.p(2,:))',1e-3)
