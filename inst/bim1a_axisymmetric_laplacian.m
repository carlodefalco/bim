## Copyright (C) 2006-2014  Carlo de Falco
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
##  author: Matteo Porro       <meoo85 _AT_ users.sourceforge.net>
##  author: Emanuela Abbate    <emanuela.abbate _AT_ mail.polimi.it>

## -*- texinfo -*-
##
## @deftypefn {Function File} @
## {@var{A}} = bim1a_axisymmetric_laplacian (@var{mesh},@var{epsilon},@var{kappa})
##
## Build the standard finite element stiffness matrix for a diffusion
## problem in cylindrical coordinates with axisymmetric configuration. 
## Rotational symmetry is assumed with respect to be the vertical
## axis r=0. Only grids that DO NOT contain r=0 are admissible.
##
##@example
##@group
##   |   |-------|   OK       |--|-----|   NO!
##  r=0                         r=0
#@end group 
#@end example
##
## The equation taken into account is:
##
## - 1/r * (r * @var{epsilon} * @var{kappa} ( u' ))' = f
## 
## where @var{epsilon} is an element-wise constant scalar function,
## while @var{kappa} is a piecewise linear conforming scalar function.
##
## @seealso{bim1a_axisymmetric_rhs, bim1a_axisymmetric_reaction, 
## bim1a_axisymmetric_advection_diffusion, bim2a_laplacian, bim3a_laplacian}
## @end deftypefn

function [A] = bim1a_axisymmetric_laplacian(mesh,epsilon,kappa)
  
  ## Check input
  if nargin != 3
    error("bim1a_axisymmetric_laplacian: wrong number of input parameters.");
  elseif !isvector(mesh)
    error("bim1a_axisymmetric_laplacian: first argument is not a valid vector.");
  endif

  ## Input-type check inside bim1a_axisymmetric_advection_diffusion
  nnodes = length(mesh);
  nelem  = nnodes - 1;
  
  A = bim1a_axisymmetric_advection_diffusion (mesh, epsilon, kappa, ones(nnodes,1), 0);
  
endfunction
