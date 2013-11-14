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
## bim2c_norm(@var{mesh},@var{u},@var{norm_type}) 
##
## Compute the @var{norm_type}-norm of function @var{u} on the domain described
## by the triangular grid @var{mesh}.
##
## The input function @var{u} can be either a piecewise linear conforming scalar
## function or an elementwise constant scalar or vector function.
##
## The string parameter @var{norm_type} can be one among 'L2', 'H1' and 'inf'.
##
## Should the input function be piecewise constant, the H1 norm will not be 
## computed and the function will return an error message.
##
## For the numerical integration of the L2 norm the second order middle point
## quadrature rule is used.
##
## @seealso{bim1c_norm, bim3c_norm}
##
## @end deftypefn

function [norm_u] = bim2c_norm (m, u, norm_type)

  ## Check input  
  if (nargin != 3)
    error ("bim2c_norm: wrong number of input parameters.");
  elseif (! (isstruct (m) && isfield (m,"p")
          && isfield (m, "t") 
          && isfield (m, "e")))
    error ("bim2c_norm: first input is not a valid mesh structure.");
  endif

  nnodes = columns (m.p);
  nel    = columns (m.t);
  
  if (isequal (size (u), [2, nel]))
    u = u';
  endif
  
  if ((length (u) != nnodes) && (rows (u) != nel))
    error ("bim2c_norm: length(u) != nnodes and rows(u) != nel.");
  endif
  
  if (! (strcmp (norm_type,'L2') 
         || strcmp (norm_type,'inf') 
         || strcmp (norm_type,'H1'))) 
    error ("bim2c_norm: invalid norm type parameter.");
  endif

  if (strcmp (norm_type,'inf'))  
    norm_u = max (abs (u(:)));
  else
    if (length (u) == nnodes)

      M = __mass_matrix__ (m);
      
      if (strcmp (norm_type, 'H1'))
        A = bim2a_laplacian (m, 1, 1);
        M += A;
      endif

      norm_u = sqrt(u' * M * u);
    
    else

      if (strcmp (norm_type, 'H1'))
        error (["bim2c_norm: cannot compute the H1 norm ", ... 
                "of an elementwise constant function."]);
      endif
      
      norm_u = m.area' * (norm (u, 2, 'rows').^2);      
      norm_u = sqrt (norm_u);

    endif
  endif

endfunction

function M = __mass_matrix__ (mesh)

  t      = mesh.t;
  nnodes = columns (mesh.p);
  nelem  = columns (t);

  ## Local contributions
  
  Mref = 1/12 * [2 1 1; 1 2 1; 1 1 2];
  area = reshape (mesh.area, 1, 1, nelem);
  
  ## Computation
  for inode = 1:3
    for jnode = 1:3
      ginode(inode,jnode,:) = t(inode,:);
      gjnode(inode,jnode,:) = t(jnode,:);
    endfor
  endfor
  Mloc = area .* Mref;

  ## assemble global matrix
  M = sparse (ginode(:), gjnode(:), Mloc(:), nnodes, nnodes);

endfunction

%!test
%!shared L, V, x, y, m
%! L = rand (1); V = rand (1); x = linspace (0,L,4);  y = x;
%! m = msh2m_structured_mesh (x,y,1,1:4);
%! m.area = msh2m_geometrical_properties (m, 'area');
%! m.shg  = msh2m_geometrical_properties (m, 'shg');
%! u    = V * ones (columns(m.p),1);
%! uinf = bim2c_norm (m, u, 'inf');
%! uL2  = bim2c_norm (m, u, 'L2');
%! uH1  = bim2c_norm (m, u, 'H1');
%! assert ([uinf, uL2, uH1], [V, V*L, V*L], 1e-12);
%!test
%! u    = V * (m.p(1,:) + 2*m.p(2,:))';
%! uinf = bim2c_norm (m, u, 'inf');
%! uL2  = bim2c_norm (m, u, 'L2');
%! uH1  = bim2c_norm (m, u, 'H1');
%! assert ([uinf, uL2, uH1], 
%!         [3*L*V, V*L^2*sqrt(8/3), V*sqrt(8/3*L^4 + 5*L^2)],
%!          1e-12);
%!test
%! u    = V * ones (columns(m.t),1);
%! uinf = bim2c_norm (m, u, 'inf');
%! uL2  = bim2c_norm (m, u, 'L2');
%! assert ([uinf, uL2], [V, V*L], 1e-12);
%!test
%! u     = V * ones (columns(m.t),1);
%! uvect = [u, 2*u];
%! uinf  = bim2c_norm (m, uvect, 'inf');
%! uL2   = bim2c_norm (m, uvect, 'L2');
%! assert ([uinf, uL2], [2*V, V*L*sqrt(5)], 1e-12);
