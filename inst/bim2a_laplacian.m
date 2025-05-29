## Copyright (C) 2006,2007,2008,2009,2010  Carlo de Falco, Massimiliano Culpo,  Copyright (C) 2025 Carlo de Falco
## 
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
##  author: Carlo de Falco     
##  author: Massimiliano Culpo 

## -*- texinfo -*-
##
## @deftypefn {Function File} @
## {@var{A}} = bim2a_laplacian (@var{mesh},@var{epsilon},@var{kappa})
##
## Build the standard finite element stiffness matrix for a diffusion
## problem. 
##
## The equation taken into account is:
##
## - div (@var{epsilon} * @var{kappa} grad (u)) = f
## 
## where @var{epsilon} is an element-wise constant scalar function,
## while @var{kappa} is a piecewise linear conforming scalar function.
##
## @seealso{bim2a_rhs, bim2a_reaction, bim2a_advection_diffusion, bim1a_laplacian, bim3a_laplacian}
## @end deftypefn

function A = bim2a_laplacian (mesh, epsilon, kappa)

  ## Check input
  if nargin != 3
    error("bim2a_laplacian: wrong number of input parameters.");
  elseif !(isstruct(mesh)     && isfield(mesh,"p") &&
	   isfield (mesh,"t") && isfield(mesh,"e"))
    error("bim2a_laplacian: first input is not a valid mesh structure.");
  endif

  p      = mesh.p;
  t      = mesh.t;
  nnodes = columns(p);
  nelem  = columns(t);

  ## Turn scalar input to a vector of appropriate size
  if isscalar(epsilon)
    epsilon  = epsilon * ones (nelem, 1);
  endif
  if isscalar(kappa)
    kappa = kappa * ones (nnodes, 1);
  endif

  if !( isvector (epsilon) && isvector (kappa) && (numel (epsilon) == nelem) && (numel (kappa) == nnodes))
    error("bim2a_laplacian: coefficients are vectors of correct size.");
  endif

  ## Local element matrices (one 3x3 matrix for each triangle in the mesh)
  Lloc = zeros(3, 3, nelem);

  ## To integrate constants over each triangle we need to multiply by
  ## the triangle area. Multiply by the diffusion coefficient now tu
  ## simplify subsequent computations. Use inverse average for the
  ## diffusion coefficient.
  kappaepsilonareak = reshape (3./sum (1./kappa(:)(mesh.t (1:3)), 1)(:) .* epsilon(:) .* mesh.area(:), 1, 1, nelem);
  shg = mesh.shg(:,:,:);
  
  ## Computation
  for inode = 1:3
    for jnode = 1:3
      ginode(inode,jnode,:) = mesh.t(inode,:);
      gjnode(inode,jnode,:) = mesh.t(jnode,:);
      Lloc(inode,jnode,:)   = sum (shg(:,inode,:) .* shg(:,jnode,:), 1) .* kappaepsilonareak;
    endfor
  endfor

  ## Assemble the local (full) matrices into one global (sparse) matrix
  A = sparse (ginode(:), gjnode(:), Lloc(:));

endfunction

%!test
%! m = msh2m_structured_mesh (0:.1:5, 0:.1:1, 1, 1:4, 'random');
%! m = bim2c_mesh_properties (m);
%! A = bim2a_laplacian (m, 1, 3);
%! dnodes = bim2c_unknowns_on_side (m, 1:4);
%! inodes = setdiff (1:columns(m.p), dnodes);
%! u = m.p(1,:)';
%! u(inodes) = A(inodes, inodes) \ (-A(inodes, dnodes) * u(dnodes));
%! assert (u, m.p(1,:)', sqrt(eps))


%!demo
%! m = msh2m_structured_mesh (0:.1:2*pi, 0:.1:2*pi, 1, 1:4, 'random');
%! m = bim2c_mesh_properties (m);
%! kappa = 2 + sin (m.p(1, :)');
%! f     = kappa .* cos (m.p(2, :)');
%! uex   = cos (m.p(2, :)');
%! A = bim2a_laplacian (m, 3./(sum(1./kappa(m.t(1:3, :)), 1)), 1);
%! b = bim2a_rhs (m, 1, f);
%! dnodes = bim2c_unknowns_on_side (m, 1:4);
%! inodes = setdiff (1:columns(m.p), dnodes);
%! u = uex;
%! u(inodes) = A(inodes, inodes) \ (b(inodes)-A(inodes, dnodes) * u(dnodes));
%! h = pdesurf (m.p, m.t, u)
%! figure
%! pdesurf (m.p, m.t, uex)
