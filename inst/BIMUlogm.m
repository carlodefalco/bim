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
## @deftypefn {Function File} @
## {[@var{T}]} = BIMUlogm (@var{t1},@var{t2})
## 
## Input:
## @itemize @minus
## @item @var{t1}:
## @item @var{t2}:
## @end itemize
##
## Output:
## @itemize @minus
## @item @var{T}:
## @end itemize
##
## @seealso{BIMUbern}
## @end deftypefn

function [T] = BIMUlogm(t1,t2)

  T = zeros(size(t2));
  
  sing     = abs(t2-t1)< 100*eps ;
  T(sing)  = (t2(sing)+t1(sing))/2;
  T(~sing) = (t2(~sing)-t1(~sing))./log(t2(~sing)./t1(~sing));

endfunction
