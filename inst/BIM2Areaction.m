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
## @deftypefn {Function File} {[@var{C}]} = BIM2Areaction(@var{mesh}, @var{delta}, @var{zeta})
##
## Builds the matrix for the discretization of the LHS
## of the equation:
## @iftex 
## @tex
## $ \delta \zeta u = f $
## @end tex 
## @end iftex 
## @ifinfo
## @var{delta} * @var{zeta} * u = f
## @end ifinfo
## 
## Input:
## @itemize @minus
## @item @var{mesh}: PDEtool-like mesh structure with required fields "p", "e", "t".
## @item @var{delta}: element-wise constant scalar function.
## @item @var{zeta}: piecewise linear conforming scalar function.
## @end itemize 
##
## @seealso{BIM2Arhs, BIM2Aadvdiff, BIM2Cmeshproperties}
## @end deftypefn

function [C] = BIM2Areaction(mesh,delta,zeta)

  Nnodes    = size(mesh.p,2);
  Nelements = size(mesh.t,2);
  
  wjacdet   = mesh.wjacdet(:,:);
  coeff     = zeta(mesh.t(1:3,:));
  coeffe    = delta(:);
  
  ## Local matrix	
  Blocmat = zeros(3,Nelements);	
  for inode = 1:3
    Blocmat(inode,:) = coeffe'.*coeff(inode,:).*wjacdet(inode,:);
  endfor
  
  gnode = (mesh.t(1:3,:));
  
  ## Global matrix
  C = sparse(gnode(:),gnode(:),Blocmat(:));

endfunction
