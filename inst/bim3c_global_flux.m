## Copyright (C) 2012  Carlo de Falco
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

## -*- texinfo -*-
##
## @deftypefn {Function File} @
## {[@var{F}]} = @
## bim3c_global_flux (@var{mesh}, @var{u}, @var{alpha}, @var{v})
##
## Compute the flux associated with the Scharfetter-Gummel approximation
## of the scalar field @var{u}.
##
## The vector field is defined as:
##
## F =- @var{alpha} ( grad (u) - grad (@var{v}) u ) 
##
## where @var{v} is a piecewise linear continuous scalar
## functions and @var{alpha} is a piecewise constant scalar function.
##
## @seealso{bim3a_rhs, bim3a_reaction, bim3a_laplacian, bim3c_mesh_properties}
## @end deftypefn


function F = bim3c_global_flux (mesh, u, acoeff, v)

  t     = mesh.t;
  nelem = columns (mesh.t);
  F     = zeros (3, nelem);

  ## Local contributions
  Lloc = zeros (4,4,nelem);

  epsilonareak = reshape (acoeff .* mesh.area', 1, 1, nelem);
  shg = mesh.shg(:,:,:);

  ## Computation
  for inode = 1:4
    for jnode = 1:4
      ginode(inode,jnode,:) = t(inode,:);
      gjnode(inode,jnode,:) = t(jnode,:);
      Lloc(inode,jnode,:)   = sum (shg(:,inode,:) .* shg(:,jnode,:), 1) ...
                              .* epsilonareak;
    endfor
  endfor

  uloc = u(t(1:4, :));
  vloc = v(t(1:4, :));
  [bp12,bm12] = bimu_bernoulli (vloc(2,:)-vloc(1,:));
  [bp13,bm13] = bimu_bernoulli (vloc(3,:)-vloc(1,:));
  [bp14,bm14] = bimu_bernoulli (vloc(4,:)-vloc(1,:));
  [bp23,bm23] = bimu_bernoulli (vloc(3,:)-vloc(2,:));
  [bp24,bm24] = bimu_bernoulli (vloc(4,:)-vloc(2,:));
  [bp34,bm34] = bimu_bernoulli (vloc(4,:)-vloc(3,:));
  
  bp12 = reshape (bp12, 1, 1, nelem) .* Lloc(1,2,:);
  bm12 = reshape (bm12, 1, 1, nelem) .* Lloc(1,2,:);
  bp13 = reshape (bp13, 1, 1, nelem) .* Lloc(1,3,:);
  bm13 = reshape (bm13, 1, 1, nelem) .* Lloc(1,3,:);
  bp14 = reshape (bp14, 1, 1, nelem) .* Lloc(1,4,:);
  bm14 = reshape (bm14, 1, 1, nelem) .* Lloc(1,4,:);
  bp23 = reshape (bp23, 1, 1, nelem) .* Lloc(2,3,:);
  bm23 = reshape (bm23, 1, 1, nelem) .* Lloc(2,3,:);
  bp24 = reshape (bp24, 1, 1, nelem) .* Lloc(2,4,:);
  bm24 = reshape (bm24, 1, 1, nelem) .* Lloc(2,4,:);
  bp34 = reshape (bp34, 1, 1, nelem) .* Lloc(3,4,:);
  bm34 = reshape (bm34, 1, 1, nelem) .* Lloc(3,4,:);
  
  ## SGloc=[...
  ##        -bm12-bm13-bm14,bp12            ,bp13           ,bp14     
  ##        bm12           ,-bp12-bm23-bm24 ,bp23           ,bp24
  ##        bm13           ,bm23            ,-bp13-bp23-bm34,bp34
  ##        bm14           ,bm24            ,bm34           ,-bp14-bp24-bp34
  ##        ];
  
  Sloc(1,1,:) = -bm12-bm13-bm14;
  Sloc(1,2,:) = bp12;
  Sloc(1,3,:) = bp13;
  Sloc(1,4,:) = bp14;

  Sloc(2,1,:) = bm12;
  Sloc(2,2,:) = -bp12-bm23-bm24; 
  Sloc(2,3,:) = bp23;
  Sloc(2,4,:) = bp24;

  Sloc(3,1,:) = bm13;
  Sloc(3,2,:) = bm23;
  Sloc(3,3,:) = -bp13-bp23-bm34;
  Sloc(3,4,:) = bp34;
  
  Sloc(4,1,:) = bm14;
  Sloc(4,2,:) = bm24;
  Sloc(4,3,:) = bm34;
  Sloc(4,4,:) = -bp14-bp24-bp34;

  r = zeros (4, nelem);
  f = zeros (3, nelem);
  
  for iel = 1:nelem

   r(:,iel) = Sloc(:,:,iel) * uloc(:,iel);
   f(:,iel) = Lloc(1:3, 1:3, iel) \ r(1:3, iel);

   F(:,iel) = shg(:,1:3, iel) * f(:, iel);

  endfor

endfunction

%!test
%! N = 10; pp = linspace (0, 1, N); msh = bim3c_mesh_properties (msh3m_structured_mesh (pp, pp, pp, 1, 1:6));
%! u = ones (N^3, 1);
%! v = ones (N^3, 1);
%! alpha = ones (columns (msh.t), 1);
%! F =  bim3c_global_flux (msh, u, alpha, v);
%! assert (norm (F(:), inf), 0, 100*eps);
