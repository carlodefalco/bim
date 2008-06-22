## Copyright (C) 2007,2008  Carlo de Falco, Massimiliano Culpo
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
##  Carlo de Falco
##  Dublin City University
##  Glasnevin, Dublin 9, Ireland
##
##  Culpo Massimiliano
##  Bergische Universitaet Wuppertal
##  Fachbereich C - Mathematik und Naturwissenschaften
##  Arbeitsgruppe fuer Angewandte MathematD-42119 Wuppertal  Gaussstr. 20 
##  D-42119 Wuppertal, Germany

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{b}]} = BIM2Arhs(@var{mesh}, @var{f}, @var{g})
##
## Constructs the RHS for the DAR problem:
## @iftex 
## @tex
##  $ -( \varepsilon  \gamma  ( u' ))' = f g$
## @end tex 
## @end iftex 
## @ifinfo
## - ( epsilon ( u' ))' =  @var{f}*@var{g}
## @end ifinfo
## 
## Input:
## @itemize @minus
## @item @var{mesh}: list of mesh nodes coordinates
## @item @var{g}: elemental values of a piecewise-wise constant function.
## @item @var{f}: nodal values of a piecewise linear conforming function.
## @end itemize 
##
## @seealso{BIM1Areaction, BIM1Alaplacian}
## @end deftypefn

function b = BIM1Arhs(mesh,f,g)

  h = (mesh(2:end)-mesh(1:end-1)).*g;
  b = f.*[h(1)/2; (h(1:end-1)+h(2:end))/2; h(end)/2];

endfunction