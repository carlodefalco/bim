## Copyright (C) 2006-2013  Carlo de Falco
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
## bim1c_norm(@var{mesh},@var{u},@var{norm_type}) 
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
## @seealso{bim2c_norm, bim3c_norm}

## @end deftypefn

function [norm_u] = bim1c_norm (m, u, norm_type)

  ## Check input  
  if (nargin != 3)
    error ("bim1c_norm: wrong number of input parameters.");
  elseif (! isvector (m))
    error ("bim1c_norm: first input is not a valid mesh.");
  endif

  nnodes = numel (m);
  nel    = numel (m) - 1;
  
  if (isrow (u))
    u = u';
  endif
  if (isrow (m))
    m = m';
  endif

  if ((numel (u) != nnodes) && (numel (u) != nel))
    error ("bim1c_norm: numel(u) != nnodes and numel(u) != nel.");
  endif
  
  if (! (strcmp (norm_type, 'L2') 
         || strcmp (norm_type, 'inf') 
         || strcmp (norm_type, 'H1')))
    error ("bim1c_norm: invalid norm type parameter.");
  endif

  if (strcmp (norm_type,'inf'))  
    norm_u = max (abs (u));
  else
    if (length (u) == nnodes)

      M = __mass_matrix__ (m);
      
      if (strcmp (norm_type, 'H1'))
        A = bim1a_laplacian (m, 1, 1);
        M += A;
      endif

      norm_u = sqrt(u' * M * u);
    
    else

      if (strcmp (norm_type, 'H1'))
        error (["bim1c_norm: cannot compute the ", ...
                "H1 norm of an elementwise constant function."]);
      endif
      
      norm_u = diff(m)' * u.^2;      
      norm_u = sqrt (norm_u);

    endif
  endif

endfunction

function M = __mass_matrix__ (m)
  
  nnodes = numel(m);
  
  h   = diff(m);
  d0  = 1/3*[h(1); h(1:end-1)+h(2:end); h(end)];
  d1  = [0; 1/6*h];
  dm1 = [1/6*h; 0];

  M  = spdiags([dm1 d0 d1], -1:1, nnodes, nnodes);
  
endfunction

%!test
%!shared L, V, m
%! L = rand (1); V = rand (1); m = linspace (0,1,5).^2;  m *= L;
%! u    = V * ones (size (m))';
%! uinf = bim1c_norm (m, u, 'inf');
%! uL2  = bim1c_norm (m, u, 'L2');
%! uH1  = bim1c_norm (m, u, 'H1');
%! assert ([uinf, uL2, uH1], [V, V*sqrt(L), V*sqrt(L)], 1e-12);
%!test
%! u    = V * m';
%! uinf = bim1c_norm (m, u, 'inf');
%! uL2  = bim1c_norm (m, u, 'L2');
%! uH1  = bim1c_norm (m, u, 'H1');
%! assert ([uinf, uL2, uH1], 
%!         [L*V, V*sqrt(L^3/3), V*sqrt(L^3/3 + L)],
%!          1e-12);
%!test
%! u    = V * ones (size (diff (m)))';
%! uinf = bim1c_norm (m, u, 'inf');
%! uL2  = bim1c_norm (m, u, 'L2');
%! assert ([uinf, uL2], [V, V*sqrt(L)], 1e-12);

