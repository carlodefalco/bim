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
## @deftypefn {Function File} {[@var{nodelist}]} = BIM3Cunknownsonside(@var{mesh},@var{sidelist})
##
## Returns the list of the mesh nodes that lie on the specified geometrical sides.
##
## Input:
## @itemize @minus
## @item @var{mesh}: PDEtool-like mesh with required field "p", "e", "t".
## @item @var{sidelist}: list of the sides of the geometrical border.
## @end itemize 
##
## Output:
## @itemize @minus
## @item @var{nodelist}: list of the nodes that lie on the specified sides.
## @end itemize
##
## @end deftypefn
  
function [nodelist] = BIM3Cunknownsonside(mesh, sidelist)
	
  [nodelist] = MSH3Mnodesonfaces(mesh,sidelist);

endfunction

%!shared mesh
% x = y = z = linspace(0,1,2);
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6);
%!test
% assert( BIM3Cunknownsonside(mesh, 1),[1 2 5 6] )
%!test
% assert( BIM3Cunknownsonside(mesh, 2),[3 4 7 8] )
%!test
% assert( BIM3Cunknownsonside(mesh, [1 2]),1:8)