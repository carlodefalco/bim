## Copyright (C) 2007,2008,2009  Carlo de Falco, Massimiliano Culpo
##
##                   BIM - Box Integration Method Package for Octave
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
##
##  AUTHORS:
##
##  Carlo de Falco <cdf _AT_ users.sourceforge.net>
##
##  Culpo Massimiliano
##  Bergische Universitaet Wuppertal
##  Fachbereich C - Mathematik und Naturwissenschaften
##  Arbeitsgruppe fuer Angewandte MathematD-42119 Wuppertal  Gaussstr. 20 
##  D-42119 Wuppertal, Germany

## -*- texinfo -*-
##
## @deftypefn {Function File} {[@var{gx},@var{gy}]} = BIM2Cpdegrad(@var{mesh},@var{u})
##
## Build the P1 approximation to the gradient of a computed solution.
##
## Input:
## @itemize @minus
## @item @var{mesh}: PDEtool-like mesh with required field "p", "e", "t".
## @item @var{u}: piecewise linear conforming scalar function.
## @end itemize 
##
## @seealso{BIM2Cglobalflux}
## @end deftypefn

function [gx, gy] = BIM2Cpdegrad(mesh,u)

  shgx = reshape(mesh.shg(1,:,:),3,[]);
  gx = sum(shgx.*u(mesh.t(1:3,:)),1);
  shgy = reshape(mesh.shg(2,:,:),3,[]);
  gy = sum(shgy.*u(mesh.t(1:3,:)),1);

endfunction
