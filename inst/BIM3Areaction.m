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
## {@var{C}} = BIM3Areaction (@var{mesh}, @var{delta}, @var{zeta})
##
## Builds the mass matrix for the discretization of the LHS of the
## equation:
## 
## @iftex 
## @tex
## $ \delta \zeta u = f $
## @end tex 
## @end iftex 
## @ifinfo
## @var{delta} * @var{zeta} * @var{u} = f
## @end ifinfo
## 
## Input:
## @itemize @minus
## @item @var{mesh}: PDEtool-like mesh structure with required fields "p", "e", "t".
## @item @var{delta}: element-wise constant scalar function.
## @item @var{zeta}: piecewise linear conforming scalar function.
## @end itemize 
##
## @seealso{BIM3Arhs, BIM3Alaplacian, BIM3Cmeshproperties}
## @end deftypefn

function [C] = BIM3Areaction (mesh,delta,zeta);

  p         = mesh.p;
  t         = mesh.t;
  nnodes    = length(p);
  nelem     = length(t);
  shp       = mesh.shp;

  C = sparse(nnodes,nnodes);

  Cloc    = zeros(4,nelem);
  coeff   = zeta(mesh.t(1:4,:));
  coeffe  = delta;
  wjacdet = mesh.wjacdet;

  for inode = 1:4
    Cloc(inode,:) = coeffe'.*coeff(inode,:).*wjacdet(inode,:);
  endfor

  gnode = (mesh.t(1:4,:));

  ## Global matrix
  C = sparse(gnode(:),gnode(:),Cloc(:));

endfunction

%!shared mesh,alpha,eta,nnodes,nelem
% x = y = z = linspace(0,1,4);
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6);
% [mesh] = BIM3Cmeshproperties(mesh);
% nnodes = columns(mesh.p);
% nelem  = columns(mesh.t);
% delta  = ones(columns(mesh.t),1);
% zeta   = ones(columns(mesh.p),1);
%!test
% [C] = BIM3Areaction(mesh,delta,zeta);
% assert(size(C),[nnodes, nnodes]);
%!test
% [C1] = BIM3Areaction(mesh,3*delta,zeta);
% [C2] = BIM3Areaction(mesh,delta,3*zeta);
% assert(C1,C2);
