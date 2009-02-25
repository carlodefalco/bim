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
## @deftypefn {Function File} {[@var{b}]} = BIM3Arhs(@var{mesh}, @var{f}, @var{g})
##
## Constructs the RHS for the DAR problem:
## @iftex 
## @tex
## $ -div ( \alpha  \gamma  ( \eta \vect{\nabla} u - \vect{beta} u )) + \delta \zeta u = f g $
## @end tex 
## @end iftex 
## @ifinfo
## -div (@var{alpha} * @var{gamma} (@var{eta} grad u - @var{beta} u )) + @var{delta} * @var{zeta} u = @var{f}*@var{g}
## @end ifinfo
## 
## Input:
## @itemize @minus
## @item @var{mesh}: PDEtool-like mesh with required field "p", "e", "t".
## @item @var{g}: element-wise constant scalar function.
## @item @var{f}: piecewise linear conforming scalar function.
## @end itemize 
##
## @seealso{BIM3Areaction, BIM3Alaplacian, BIM3Cmeshproperties}
## @end deftypefn

function [b] = BIM3Arhs (mesh,f,g);

  p         = mesh.p;
  t         = mesh.t;
  nnodes    = length(p);
  nelem     = length(t);
  shp       = mesh.shp;

  b = sparse(nnodes,1);

  bloc    = zeros(4,nelem);
  coeff   = f(mesh.t(1:4,:));
  coeffe  = g;
  wjacdet = mesh.wjacdet;

  for inode = 1:4
    bloc(inode,:) = coeffe'.*coeff(inode,:).*wjacdet(inode,:);
  endfor

  gnode = (mesh.t(1:4,:));

  ## Global matrix
  b = sparse(gnode(:),1,bloc(:));

endfunction

%!shared mesh,alpha,eta,nnodes,nelem
% x = y = z = linspace(0,1,4);
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6);
% [mesh] = BIM3Cmeshproperties(mesh);
% nnodes = columns(mesh.p);
% nelem  = columns(mesh.t);
% g      = ones(columns(mesh.t),1);
% f      = ones(columns(mesh.p),1);
%!test
% [b] = BIM3Arhs(mesh,f,g);
% assert(size(b),[nnodes, 1]);
%!test
% [b1] = BIM3Arhs(mesh,3*f,g);
% [b2] = BIM3Arhs(mesh,f,3*g);
% assert(b1,b2);
