## Copyright (C) 2013 Carlo de Falco
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
##
## @deftypefn {Function File} {@var{u_nod}} = bim1c_elem_to_nodes (@var{mesh}, @var{u_el}) 
## @deftypefnx {Function File} {@var{u_nod}} = bim1c_elem_to_nodes (@var{m_el}, @var{u_el}) 
## @deftypefnx {Function File} {[@var{u_nod}, @var{m_el}]} = bim1c_elem_to_nodes ( ... ) 
##
## Compute interpolated values at nodes @var{u_nod} given values at element mid-points @var{u_el}.
## If called with more than one output, also return the interpolation matrix @var{m_el} such that
## @code{u_nod = m_el * u_el}.
## If repeatedly performing interpolation on the same mesh the matrix @var{m_el} obtained by a previous call 
## to @code{bim1c_elem_to_nodes} may be passed as input to avoid unnecessary computations.
##
## @end deftypefn


## Author: Carlo de Falco <cdf _AT_ users.sourceforge.net>
## Author: Matteo Porro <meoo85 _AT_ users.sourceforge.net>
## Created: 2013-11-04

function [u_nod, m_el] = bim1c_elem_to_nodes (m, u_el)

  if (nargout > 1 )
    if (isvector (m))
      nel  = numel (m) - 1;
      nnod = numel (m);
      m_el = spalloc (nnod, nel, 2 * nel);
      h = diff (m);
      for iel = 1:nel
        m_el([iel, iel+1], iel) = h(iel);
      endfor
      m_el = diag (sum (m_el, 2)) \ m_el;
    elseif (ismatrix (m))
      m_el = m;
    else
      error (["bim1c_elem_to_nodes: first input ", ...
              "parameter is of incorrect type"]);
    endif
    u_nod = m_el * u_el;
  else
    if (isvector (m))
      rhs  = bim1a_rhs (m, u_el, 1);
      mass = bim1a_reaction (m, 1, 1);
      u_nod = full (mass \ rhs);
    elseif (ismatrix (m))
      u_nod = m * u_el;
    else
      error (["bim1c_elem_to_nodes: first input ", ...
              "parameter is of incorrect type"]);
    endif      
  endif

endfunction

%!test
%! n = 10; msh = linspace (0, 1, n+1);
%! nel  = n;
%! nnod = n+1;
%! u_el = randn (nel, 1);
%! un1 = bim1c_elem_to_nodes (msh, u_el);
%! [un2, m] = bim1c_elem_to_nodes (msh, u_el);
%! un3 = bim1c_elem_to_nodes (m, u_el);
%! [un4, m] = bim1c_elem_to_nodes (m, u_el);
%! assert (un1, un2, 1e-10)
%! assert (un1, un3, 1e-10)
%! assert (un1, un4, 1e-10)

