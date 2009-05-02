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
## {@var{A}} = BIM3Alaplacian (@var{mesh}, @var{alpha}, @var{eta})
##
## Builds the Laplacian matrix for the discretization of the LHS of the
## equation:
## @iftex 
## @tex
## $ -div ( \alpha  ( \eta \vect{\nabla} u ) $
## @end tex 
## @end iftex 
## @ifinfo
## -div (@var{alpha} (@var{eta} grad u ) = f
## @end ifinfo
## 
## where: 
## @itemize @minus
## @item @var{alpha}: element-wise constant scalar function
## @item @var{eta}: piecewise linear conforming scalar functions
## @end itemize
##
## @seealso{BIM3Arhs, BIM3Areaction, BIM3Cmeshproperties}
## @end deftypefn

function [A] = BIM3Alaplacian (mesh,alpha,eta)

  p      = mesh.p;
  t      = mesh.t;
  nnodes = columns(p);
  nelem  = columns(t);

  ## Check input
  if ( min(size(alpha)) + min(size(eta)) )  != 2
    ## check if alpha or eta are matrices
    ## multidimensional array are not taken into account
    warning("Check size of alpha and beta!");
    print_usage;
  elseif length(alpha) != nelem
    warning("Length of alpha is not equal to the number of elements!");
    print_usage;
  elseif length(eta) != nnodes
    warning("Length of eta is not equal to\nthe number of nodes!");
    print_usage;
  endif

  Lloc = zeros(4,4,nelem);

  alphaareak = reshape (alpha .* mesh.area',1,1,nelem);
  shg        = mesh.shg(:,:,:);
  
  for inode = 1:4
    for jnode = 1:4
      ginode(inode,jnode,:) = mesh.t(inode,:);
      gjnode(inode,jnode,:) = mesh.t(jnode,:);
      Lloc(inode,jnode,:)   = \
	  sum( eta(inode) * shg(:,inode,:) .* shg(:,jnode,:),1) .* alphaareak;
    endfor
  endfor

  A = sparse(ginode(:),gjnode(:),Lloc(:));

endfunction

%!shared mesh,alpha,eta,nnodes,nelem
% x = y = z = linspace(0,1,4);
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6);
% [mesh] = BIM3Cmeshproperties(mesh);
% nnodes = columns(mesh.p);
% nelem  = columns(mesh.t);
% alpha = ones(columns(mesh.t),1);
% eta   = ones(columns(mesh.p),1);
%!test
% [A] = BIM3Alaplacian(mesh,alpha,eta);
% assert(size(A),[nnodes, nnodes]);
%!test
% [A1] = BIM3Alaplacian(mesh,3*alpha,eta);
% [A2] = BIM3Alaplacian(mesh,alpha,3*eta);
% assert(A1,A2);