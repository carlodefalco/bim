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
## @deftypefn {Function File} {[@var{omesh}]} = BIM3Cmeshproperties(@var{imesh})
##
## Compute mesh properties needed by the operator assembly routines and
## add them as fields to the mesh structure.
##
## Input:
## @itemize @minus
## @item @var{imesh}: PDEtool-like mesh with required field "p", "e", "t".
## @end itemize 
##
## Output:
## @itemize @minus
## @item @var{omesh}: PDEtool-like mesh structure with added fields needed by BIM method.
## @end itemize
##
## @seealso{BIM3Areaction, BIM3Alaplacian, BIM3Arhs}
## @end deftypefn

function [omesh] = BIM3Cmeshproperties(imesh)

  omesh = imesh;
  [omesh.wjacdet,omesh.area,omesh.shg,omesh.shp] = \
      MSH3Mgeomprop(imesh,"wjacdet","area","shg","shp");

endfunction

%!shared mesh
% x = y = z = linspace(0,1,4);
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6);
% [mesh] = BIM3Cmeshproperties(mesh);
%!test
% tmp = MSH3Mgeomprop(mesh,"wjacdet");
% assert(mesh.wjacdet,tmp);
%!test
% tmp = MSH3Mgeomprop(mesh,"shg");
% assert(mesh.shg,tmp);
%!test
% tmp = MSH3Mgeomprop(mesh,"shp");
% assert(mesh.shp,tmp);
%!test
% assert(mesh.area,sum(mesh.wjacdet,1));

