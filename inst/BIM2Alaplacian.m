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
##
## @deftypefn {Function File} @
## {@var{A}} = BIM2Alaplacian (@var{mesh}, @var{epsilon})
##
## Builds the finite-element matrix for the 
## discretization of the LHS
## of the equation:
## 
## @iftex 
## @tex
## $ -div ( \varepsilon  \gamma  ( \nabla u )) = f $
## @end tex 
## @end iftex 
## @ifinfo
## - div (@var{epsilon} grad ( u )) = f
## @end ifinfo
## 
## where: 
## @itemize @minus
## @item @var{epsilon}: elemental values of an piece-wise constant function
## @end itemize
##
##
## @seealso{BIM2Arhs, BIM2Areaction}
## @end deftypefn

function A = BIM2Alaplacian(mesh,epsilon)
  Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
  A = BIM2Aadvdiff (mesh, epsilon, ones(Nnodes,1), ones(Nnodes,1), 0);
endfunction