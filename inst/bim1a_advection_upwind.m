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
##  author: Matteo Porro       <meoo85 _AT_ users.sourceforge.net>

## -*- texinfo -*-
##
## @deftypefn {Function File} @
## {[@var{A}]} = bim1a_advection_upwind (@var{mesh}, @var{beta})
##
## Build the UW stabilized stiffness matrix for an advection problem. 
## 
## The equation taken into account is:
##
##  (@var{beta} u)' = f
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
## @seealso{bim1a_rhs, bim1a_reaction, bim1a_laplacian, bim2a_advection_diffusion} 
## @end deftypefn

function A = bim1a_advection_upwind (x, beta)

  ## Check input
  if nargin != 2
    error("bim1a_advection_upwind: wrong number of input parameters.");
  endif
  
  nnodes = length(x);
  nelem  = nnodes-1;

  if (length(beta) == 1)
    vk = zeros(nelem,1);
  elseif (length(beta) == nelem)
    vk = beta; 
  elseif (length(beta) == nnodes)
    vk = diff(beta);
  else
    error("bim1a_advection_upwind: coefficient beta has wrong dimensions.");
  endif
  
  bmk =  (vk+abs(vk))/2;
  bpk = -(vk-abs(vk))/2;
 
  dm1 = [-bmk; NaN];
  dp1 = [NaN; -bpk]; 
  d0  = [bmk(1); bmk(2:end) + bpk(1:end-1); bpk(end)];
  A   = spdiags([dm1, d0, dp1],-1:1,nnodes,nnodes);

endfunction

%!test
%! n = 200;
%! mesh      = linspace(0,1,n+1)';
%! uex       = @(r) - r.^2 + 1;
%! Nnodes    = numel(mesh);
%! Nelements = Nnodes-1;
%! D = 1;  v = 1;  sigma = 0;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 1/D*v*ones(Nelements,1);
%! delta  = ones(Nelements,1);
%! zeta   = sigma*ones(Nnodes,1);
%! f      = @(r) 2*D - 2*v.*r + sigma*uex(r);
%! rhs    = bim1a_rhs(mesh, ones(Nelements,1), f(mesh));
%! S = bim1a_laplacian(mesh,alpha,gamma);
%! A = bim1a_advection_upwind(mesh, beta);
%! R = bim1a_reaction(mesh, delta, zeta);
%! S += (A+R);
%! u = zeros(Nnodes,1); u([1 end]) = uex(mesh([1 end]));
%! u(2:end-1) = S(2:end-1,2:end-1)\(rhs(2:end-1) - S(2:end-1,[1 end])*u([1 end]));
%! assert(u,uex(mesh),1e-3)
