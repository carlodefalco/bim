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
## @deftypefn {Function File} {[@var{C}]} = @
## bim1a_axisymmetric_reaction(@var{mesh},@var{delta},@var{zeta})
##
## Build the lumped finite element mass matrix for a diffusion
## problem in cylindrical coordinates with axisymmetric configuration. 
##
## The equation taken into account is:
##
## @var{delta} * @var{zeta} * u = f
## 
## where @var{delta} is an element-wise constant scalar function, while
## @var{zeta} is a piecewise linear conforming scalar function.
##
## @seealso{bim1a_axisymmetric_rhs, bim1a_axisymmetric_advection_diffusion, bim1a_axisymmetric_laplacian,
## bim2a_reaction, bim3a_reaction}
## @end deftypefn

function [C] = bim1a_axisymmetric_reaction(mesh,delta,zeta)
  
  ## Check input
  if nargin != 3
    error("bim1a_axisymmetric_reaction: wrong number of input parameters.");
  elseif !isvector(mesh)
    error("bim1a_axisymmetric_reaction: first argument is not a valid vector.");
  endif

  mesh    = reshape(mesh,[],1);
  nnodes  = length(mesh);
  nelems  = nnodes-1;

  ## Turn scalar input to a vector of appropriate size
  if isscalar(delta)
    delta = delta*ones(nelems,1);
  endif
  if isscalar(zeta)
    zeta  = zeta*ones(nnodes,1);
  endif

  if !( isvector(delta) && isvector(zeta) )
    error("bim1a_axisymmetric_reaction: coefficients are not valid vectors.");
  elseif length(delta) != nelems
    error("bim1a_axisymmetric_reaction: length of delta is not equal to the number of elements.");
  elseif length(zeta)  != nnodes
    error("bim1a_axisymmetric_reaction: length of zeta is not equal to the number of nodes.");
  endif

  h 	= (mesh(2:end)-mesh(1:end-1)).*delta;
  d0	= zeta.*[h(1)/2; (h(1:end-1)+h(2:end))/2; h(end)/2];
  C   = spdiags(d0.*abs(mesh), 0, nnodes,nnodes);

endfunction

%!test
%! n = 100;
%! mesh      = linspace(0,1,n+1)';
%! cm        = (mesh(1:end-1) + mesh(2:end))/2; 
%! uex       = @(r) - r.^2 + 1;
%! Nnodes    = numel(mesh);
%! Nelements = Nnodes-1;
%! D = 1;  v = cm;  sigma = 1;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 0;
%! delta  = ones(Nelements,1);
%! zeta   = sigma*ones(Nnodes,1);
%! f      = @(r) 4*D + sigma*uex(r);
%! rhs    = bim1a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh));
%! S = bim1a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! R = bim1a_axisymmetric_reaction(mesh, delta, zeta);
%! S += R;
%! u = zeros(Nnodes,1); u(end) = uex(mesh(end));
%! u(1:end-1) = S(1:end-1,1:end-1)\(rhs(1:end-1) - S(1:end-1,end)*u(end));
%! assert(u,uex(mesh),1e-3)

%!test
%! n = 100;
%! mesh      = linspace(0,1,n+1)';
%! cm        = (mesh(1:end-1) + mesh(2:end))/2; 
%! uex       = @(r) - r.^2 + 1;
%! Nnodes    = numel(mesh);
%! Nelements = Nnodes-1;
%! D = 1;  v = cm;  sigma = 1;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 1/D*v;
%! delta  = ones(Nelements,1);
%! zeta   = sigma*ones(Nnodes,1);
%! f      = @(r) 4*D + 2 - 4*r.^2 + sigma*uex(r);
%! rhs    = bim1a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh));
%! S = bim1a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! R = bim1a_axisymmetric_reaction(mesh, delta, zeta);
%! S += R;
%! u = zeros(Nnodes,1); u(end) = uex(mesh(end));
%! u(1:end-1) = S(1:end-1,1:end-1)\(rhs(1:end-1) - S(1:end-1,end)*u(end));
%! assert(u,uex(mesh),1e-3)

%!test
%! x = linspace(0,1,101);
%! A = bim1a_axisymmetric_reaction(x,1,1);
%! delta = ones(100,1);
%! zeta  = ones(101,1);
%! B = bim1a_axisymmetric_reaction(x,delta,zeta);
%! assert(A,B)
