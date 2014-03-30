## Copyright (C) 2006,2007,2008,2009,2010  Carlo de Falco, Massimiliano Culpo
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
##  author: Carlo de Falco <cdf _AT_ users.sourceforge.net>
##  author: Matteo Porro   <meoo85 _AT_ users.sourceforge.net>

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{M}]} = @
## bim3a_boundary_mass(@var{mesh},@var{facelist},@var{nodelist})
##
## Build the lumped boundary mass matrix needed to apply Robin boundary
## conditions.
##
## The vector @var{facelist} contains the list of the faces contributing
## to the mass matrix.
##
## The optional argument @var{nodelist} contains the list of the
## degrees of freedom on the boundary.
##
## @seealso{bim3a_rhs, bim3a_advection_diffusion, bim3a_laplacian,
## bim3a_reaction, bim2a_boundary_mass} 
## @end deftypefn

function [M] = bim3a_boundary_mass (mesh, facelist, nodelist)

  ## Check input
  if (nargin > 3)
    error ("bim3a_boundary_mass: wrong number of input parameters.");
  elseif (! ((isstruct (mesh)) && (isfield (mesh, "p")) 
             && (isfield (mesh, "t")) && isfield(mesh, "e")))
    error (["bim3a_boundary_mass: first input", ...
            " is not a valid mesh structure."]);
  elseif (! ((isvector (facelist)) && (isnumeric (facelist))))
    error (["bim3a_boundary_mass: second ", ...
            "input is not a valid numeric vector."]);
  endif

  if (nargin < 3)
    [nodelist] = bim3c_unknowns_on_faces (mesh, facelist);
  endif

  p = mesh.p;
  t = [];
  for ie = facelist
    t = [t,  mesh.e([1:3 10], mesh.e(10,:) == ie)];
  endfor

  area = 1/2 * norm (cross (p(:,t(2,:))-p(:,t(1,:)), 
                            p(:,t(3,:))-p(:,t(1,:))), 
                     2, 'columns');
  
  dd = zeros (size (nodelist));
  
  for in = 1:numel (nodelist)
    dd (in) = 1/3 * sum (area (any (t(1:3,:) == nodelist(in))));
  endfor
  
  M = sparse (diag (dd));

endfunction

