## Copyright (C) 2006-2014 Carlo de Falco, Massimiliano Culpo
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
##
## @deftypefn {Function File} @
## {[@var{A}]} = @
## bim2a_axisymmetric_advection_upwind (@var{mesh}, @var{beta})
##
## Build the Upwind stabilized stiffness matrix for an advection problem 
## in cylindrical coordinates with axisymmetric configuration.
## 
## The equation taken into account is:
##
## 1/r * d/dr (r * @var{beta}_r u) + d/dz (@var{beta}_z u) = f
##
## where @var{beta} is an element-wise constant vector function.
##
## Instead of passing the vector field @var{beta} directly one can pass
## a piecewise linear conforming scalar function  @var{phi} as the last
## input.  In such case @var{beta} = grad @var{phi} is assumed.
##
## If @var{phi} is a single scalar value @var{beta} is assumed to be 0
## in the whole domain.
## 
## @seealso{bim2a_axisymmetric_rhs, bim2a_axisymmetric_reaction, 
## bim2a_axisymmetric_advection_diffusion, bim2c_mesh_properties}
## @end deftypefn

function A = bim2a_axisymmetric_advection_upwind (mesh, beta)

  ## Check input
  if nargin != 2
    error("bim2a_axisymmetric_advection_upwind: wrong number of input parameters.");
  elseif !(isstruct(mesh)     && isfield(mesh,"p") &&
	   isfield (mesh,"t") && isfield(mesh,"e"))
    error("bim2a_axisymmetric_advection_upwind: first input is not a valid mesh structure.");
  elseif !(all(mesh.p(1,:) >= 0) || all(mesh.p(1,:) <= 0))
    error("bim2a_axisymmetric_advection_upwind: the input mesh cannot intersect the rotation axis r=0.");
  endif

  nnodes = columns(mesh.p);
  nelem  = columns(mesh.t);
  
  x = abs (mesh.p(1,:));
  x = x(mesh.t(1:3,:));
  y = mesh.p(2,:);
  y = y(mesh.t(1:3,:));
  
  alphaareak = reshape (mesh.area, 1, 1, nelem);
  shg        = mesh.shg(:,:,:);
  
  ## Build local Laplacian matrix	
  Lloc = zeros(3,3,nelem);	
  
  for inode = 1:3
    for jnode = 1:3
      ginode(inode,jnode,:) = mesh.t(inode,:);
      gjnode(inode,jnode,:) = mesh.t(jnode,:);
      Lloc(inode,jnode,:)   = sum( shg(:,inode,:) .* shg(:,jnode,:),1) .* alphaareak;
    endfor
  endfor
  
  if all(size(beta)==1)
    v12 = 0;
    v23 = 0;
    v31 = 0; 
  elseif all(size(beta)==[2,nelem])
    v12 = beta(1,:) .* (x(2,:)-x(1,:)) + beta(2,:) .* (y(2,:)-y(1,:));
    v23 = beta(1,:) .* (x(3,:)-x(2,:)) + beta(2,:) .* (y(3,:)-y(2,:));
    v31 = beta(1,:) .* (x(1,:)-x(3,:)) + beta(2,:) .* (y(1,:)-y(3,:)); 
  elseif all(size(beta)==[nnodes,1])
    betaloc = beta(mesh.t(1:3,:));
    v12     = betaloc(2,:)-betaloc(1,:);
    v23     = betaloc(3,:)-betaloc(2,:);
    v31     = betaloc(1,:)-betaloc(3,:); 
  else
    error("bim2a_axisymmetric_advection_upwind: coefficient beta has wrong dimensions.");
  endif
  
  [bp12, bm12] = deal (- (v12 - abs (v12))/2, (v12 + abs (v12))/2);
  [bp23, bm23] = deal (- (v23 - abs (v23))/2, (v23 + abs (v23))/2);
  [bp31, bm31] = deal (- (v31 - abs (v31))/2, (v31 + abs (v31))/2);
  
  r12 = (x(2,:) + x(1,:)) / 2;
  r23 = (x(3,:) + x(2,:)) / 2;
  r31 = (x(1,:) + x(3,:)) / 2;

  bp12 = reshape(r12 .* bp12,1,1,nelem).*Lloc(1,2,:);
  bm12 = reshape(r12 .* bm12,1,1,nelem).*Lloc(1,2,:);
  bp23 = reshape(r23 .* bp23,1,1,nelem).*Lloc(2,3,:);
  bm23 = reshape(r23 .* bm23,1,1,nelem).*Lloc(2,3,:);
  bp31 = reshape(r31 .* bp31,1,1,nelem).*Lloc(3,1,:);
  bm31 = reshape(r31 .* bm31,1,1,nelem).*Lloc(3,1,:);
  
  Sloc(1,1,:) = (-bm12-bp31);
  Sloc(1,2,:) = bp12;
  Sloc(1,3,:) = bm31;
  
  Sloc(2,1,:) = bm12;
  Sloc(2,2,:) = (-bp12-bm23); 
  Sloc(2,3,:) = bp23;
  
  Sloc(3,1,:) = bp31;
  Sloc(3,2,:) = bm23;
  Sloc(3,3,:) = (-bm31-bp23);
  
  A = sparse(ginode(:), gjnode(:), Sloc(:));

endfunction
