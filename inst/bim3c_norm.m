## Copyright (C) 2006-2013  Carlo de Falco, Massimiliano Culpo
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
##  author: Matteo Porro       <meoo85 _AT_ users.sourceforge.net>

## -*- texinfo -*-
##
## @deftypefn {Function File} {[@var{norm_u}]} = @
## bim3c_norm(@var{mesh},@var{u},@var{norm_type}) 
##
## Compute the @var{norm_type}-norm of function @var{u} on the domain described
## by the tetrahedral grid @var{mesh}.
##
## The input function @var{u} can be either a piecewise linear conforming scalar
## function or an elementwise constant scalar or vector function.
##
## The string parameter @var{norm_type} can be one among 'L2', 'H1' and 'inf'.
##
## Should the input function be piecewise constant, the H1 norm will not be 
## computed and the function will return an error message.
##
## For the numerical integration of the L2 norm the second order quadrature rule
## by Keast is used (ref. P. Keast, Moderate degree tetrahedral quadrature
## formulas, CMAME 55: 339-348 1986).
##
## @seealso{bim1c_norm, bim2c_norm}

## @end deftypefn

function [norm_u] = bim3c_norm (m, u, norm_type)

  ## Check input  
  if (nargin != 3)
    error ("bim3c_norm: wrong number of input parameters.");
  elseif (! (isstruct (m) && isfield (m,"p")) && isfield (m, "t") 
          && isfield (m, "e"))
    error ("bim3c_norm: first input is not a valid mesh structure.");
  endif

  nnodes = columns (m.p);
  nel    = columns (m.t);
  
  if (isequal (size (u), [3, nel]))
    u = u';
  endif
  
  if ((length (u) != nnodes) && (rows (u) != nel))
    error ("bim3c_norm: length(u) != nnodes and rows(u) != nel.");
  endif
  
  if !(strcmp (norm_type,'L2') || strcmp (norm_type,'inf') || 
       strcmp (norm_type,'H1')) 
    error ("bim3c_norm: invalid norm type parameter.");
  endif

  if (strcmp (norm_type,'inf'))  
    norm_u = max (abs (u(:)));
  else
    if (length (u) == nnodes)

      M = __mass_matrix__ (m);
      
      if (strcmp (norm_type, 'H1'))
        A = bim3a_laplacian (m, 1, 1);
        M += A;
      endif

      norm_u = sqrt(u' * M * u);
    
    else

      if (strcmp (norm_type, 'H1'))
        error ("bim3c_norm: cannot compute the H1 norm of an elementwise constant function.");
      endif
      
      norm_u = m.area * (norm (u', 2, 'cols').^2)';      
      norm_u = sqrt (norm_u);

    endif
  endif

endfunction

function M = __mass_matrix__ (mesh)

  t      = mesh.t;
  nnodes = columns (mesh.p);
  nelem  = columns (t);

  ## Local contributions
  a = (5 + 3 * sqrt (5)) / 20;   b = (5 - sqrt (5)) / 20;
  l1 = (1 - 3*b)^2 + 3*(1 - 2*b - a)^2;
  l2 = (1 - 3*b)*b + (1 - 2*b - a)*(a + 2*b);
  
  Mref = 1/4 * [l1 l2 l2 l2; l2 l1 l2 l2; l2 l2 l1 l2; l2 l2 l2 l1];
  area = reshape (mesh.area, 1, 1, nelem);
  
  ## Computation
  for inode = 1:4
    for jnode = 1:4
      ginode(inode,jnode,:) = t(inode,:);
      gjnode(inode,jnode,:) = t(jnode,:);
    endfor
  endfor	
  Mloc = area .* Mref;

  ## assemble global matrix
  M = sparse (ginode(:), gjnode(:), Mloc(:), nnodes, nnodes);

endfunction

%!test
%!shared L, V, x, y, z, m
%! L = rand (1); V = rand (1); x = linspace (0,L,4);  y = x;  z = x;
%! m = msh3m_structured_mesh (x,y,z,1,1:6);
%! m.area = msh3m_geometrical_properties (m, 'area');
%! m.shg  = msh3m_geometrical_properties (m, 'shg');
%! u    = V * ones (columns(m.p),1);
%! uinf = bim3c_norm (m, u, 'inf');
%! uL2  = bim3c_norm (m, u, 'L2');
%! uH1  = bim3c_norm (m, u, 'H1');
%! assert ([uinf, uL2, uH1], [V, V*sqrt(L^3), V*sqrt(L^3)], 1e-12);
%!test
%! u    = V * (m.p(1,:) + 2*m.p(2,:) + 3*m.p(3,:))';
%! uinf = bim3c_norm (m, u, 'inf');
%! uL2  = bim3c_norm (m, u, 'L2');
%! uH1  = bim3c_norm (m, u, 'H1');
%! assert ([uinf, uL2, uH1], 
%!         [6*L*V, V*sqrt(61/6*L^5), V*sqrt(61/6*L^5 + 14*L^3)],
%!          1e-12);
%!test
%! u    = V * ones (columns(m.t),1);
%! uinf = bim3c_norm (m, u, 'inf');
%! uL2  = bim3c_norm (m, u, 'L2');
%! assert ([uinf, uL2], [V, V*sqrt(L^3)], 1e-12);
%!test
%! u     = V * ones (columns(m.t),1);
%! uvect = [u, 2*u, 3*u];
%! uinf  = bim3c_norm (m, uvect, 'inf');
%! uL2   = bim3c_norm (m, uvect, 'L2');
%! assert ([uinf, uL2], [3*V, V*sqrt(14*L^3)], 1e-12);
