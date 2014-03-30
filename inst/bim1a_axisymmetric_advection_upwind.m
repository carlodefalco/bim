## Copyright (C) 2010-2014  Carlo de Falco
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
## {[@var{A}]} = bim1a_axisymmetric_advection_upwind (@var{mesh}, @var{beta})
##
## Build the Upwind stabilized stiffness matrix for an advection problem 
## in cylindrical coordinates with axisymmetric configuration.
## 
## The equation taken into account is:
##
## 1/r * (r * @var{beta} u)' = f
##
## where @var{beta} is an element-wise constant.
##
## Instead of passing the vector field @var{beta} directly one can pass
## a piecewise linear conforming scalar function  @var{phi} as the last
## input.  In such case @var{beta} = grad @var{phi} is assumed.
##
## If @var{phi} is a single scalar value @var{beta} is assumed to be 0
## in the whole domain. 
##
## @seealso{bim1a_axisymmetric_advection_diffusion, bim1a_axisymmetric_rhs,
## bim1a_axisymmetric_reaction, bim1a_axisymmetric_laplacian}
## @end deftypefn

function A = bim1a_axisymmetric_advection_upwind (x, beta)

  ## Check input
  if nargin != 2
    error("bim1a_axisymmetric_advection_upwind: wrong number of input parameters.");
  endif
  
  nnodes = length(x);
  nelem  = nnodes-1;

  cm    = reshape((x(1:end-1)+x(2:end))/2,[],1)
  
  if (length(beta) == 1)
    vk = 0;#zeros(nelem,1);
  elseif (length(beta) == nelem)
    vk = beta; 
  elseif (length(beta) == nnodes)
    vk = diff(beta);
  else
    error("bim1a_axisymmetric_advection_upwind: coefficient beta has wrong dimensions.");
  endif

  bmk =  (vk+abs(vk))/2 .* abs(cm);
  bpk = -(vk-abs(vk))/2 .* abs(cm);
 
  dm1 = [-bmk; NaN]; 
  dp1 = [NaN; -bpk]; 
  d0  = [bmk(1); bmk(2:end) + bpk(1:end-1); bpk(end)];
  A   = spdiags([dm1, d0, dp1],-1:1,nnodes,nnodes);

endfunction

%!test
%! nn = 20;
%! mesh = linspace(1,2,nn+1)';
%! D = 1;  v = 0;  sigma = 0;
%! uex      = @(r) exp(r);
%! duexdr   = @(r) uex(r);
%! d2uexdr2 = @(r) uex(r);
%! f = @(r,z) -D./r.*duexdr(r) - D.*d2uexdr2(r) ...
%!           + v./r .* uex(r) + v * duexdr(r) ...
%!           + sigma * uex(r);
%! uex_left = uex(mesh(1));  uex_right = uex(mesh(end));
%! Ar  = bim1a_axisymmetric_laplacian (mesh, D, 1);
%! Adv = bim1a_axisymmetric_advection_upwind (mesh, v*ones(nn,1));
%! R   = bim1a_axisymmetric_reaction (mesh, sigma, 1);
%! M = Ar + Adv + R;
%! M(1,:)   *= 0;  M(1,1) = 1;
%! M(end,:) *= 0;  M(end, end)  = 1;
%! rhs = bim1a_axisymmetric_rhs (mesh, 1, f(mesh));
%! rhs(1)   = uex_left;  rhs(end) = uex_right;
%! uh = M \ rhs;  
%! assert(uh, uex(mesh), 1e-3);
