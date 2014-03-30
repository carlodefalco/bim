## Copyright (C) 2006-2014  Carlo de Falco, Massimiliano Culpo
##
## This file is part of:
##     BIM - Diffusion Advection Reaction PDE Solver
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
##  author: Carlo de Falco     <cdf _AT_ users.sourceforge.net>
##  author: Massimiliano Culpo <culpo _AT_ users.sourceforge.net>
##  author: Matteo porro       <meoo85 _AT_ users.sourceforge.net>
##  author: Emanuela Abbate    <emanuela.abbate _AT_ mail.polimi.it>

## -*- texinfo -*-
##
## @deftypefn {Function File} @
## {@var{A}} = bim2a_axisymmetric_laplacian (@var{mesh},@var{epsilon},@var{kappa})
##
## Build the standard finite element stiffness matrix for a diffusion
## problem in cylindrical coordinates with axisymmetric configuration.
## Rotational symmetry is assumed with respect to be the vertical axis r=0. 
## Only plane geometries that DO NOT intersect the symmetry axis are admitted.
##
##@example
##@group
##    |   ____                 _|____ 
##    |  |    \               \ |    |
##  z |  |     \  OK           \|    |   NO!
##    |  |______\               |\___|
##    |     r                   |
#@end group 
#@end example
##
## The equation taken into account is:
##
## 1/r * d(r * Fr)/dr + dFz/dz = f
##
## with
##
## F = [Fr, Fz]' = - @var{epsilon} * @var{kappa} grad (u)
## 
## where @var{epsilon} is an element-wise constant scalar function,
## while @var{kappa} is a piecewise linear conforming scalar function.
##
## @seealso{bim2a_axisymmetric_rhs, bim2a_axisymmetric_reaction, 
## bim2a_axisymmetric_advection_diffusion, bim2a_laplacian, bim1a_laplacian, 
## bim3a_laplacian}
## @end deftypefn

function [A] = bim2a_axisymmetric_laplacian(mesh,epsilon,kappa)

  ## Check input
  if nargin != 3
    error("bim2a_axisymmetric_laplacian: wrong number of input parameters.");
  elseif !(all(mesh.p(1,:) >= 0) || all(mesh.p(1,:) <= 0))
    error("bim2a_axisymmetric_laplacian: the input mesh cannot intersect the rotation axis r=0.");
  endif

  ## Input check inside bim2a_axisymmetric_advection_diffusion
  nnodes = columns(mesh.p);
  nelem  = columns(mesh.t);
  
  A = bim2a_axisymmetric_advection_diffusion (mesh,epsilon,kappa,ones(nnodes,1),0);

endfunction
