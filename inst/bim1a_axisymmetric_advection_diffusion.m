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
##  author: Matteo Porro       <meoo85 _AT_ users.sourceforge.net>
##  author: Emanuela Abbate    <emanuela.abbate _AT_ mail.polimi.it>

## -*- texinfo -*-
##
## @deftypefn {Function File} @
## {[@var{A}]} = @
## bim1a_axisymmetric_advection_diffusion(@var{mesh},@var{alpha},@var{gamma},@var{eta},@var{beta})
##
## Build the Scharfetter-Gummel stabilized stiffness matrix for a
## diffusion-advection problem in cylindrical coordinates with axisymmetric 
## configuration. Rotational symmetry is assumed with respect to be the vertical
## axis r=0. Only grids that DO NOT contain r=0 are admissible.
##
##@example
##@group
##   |   |-------|   OK       |-------|   |   OK        |--|-----|   NO!
##  r=0                                  r=0              r=0
#@end group 
#@end example
## 
## The equation taken into account is:
##
## - 1/r * d/dr (@var{alpha} * @var{gamma} (@var{eta} du/dr - @var{beta} u)) = f
##
## where @var{alpha} is an element-wise constant scalar function,
## @var{eta} and @var{gamma} are piecewise linear conforming scalar
## functions, @var{beta} is an element-wise constant vector function.
##
## Instead of passing the vector field @var{beta} directly one can pass
## a piecewise linear conforming scalar function  @var{phi} as the last
## input.  In such case @var{beta} = grad @var{phi} is assumed.
##
## If @var{phi} is a single scalar value @var{beta} is assumed to be 0
## in the whole domain. 
##
## @seealso{bim1a_axisymmetric_rhs, bim1a_axisymmetric_reaction, 
## bim1a_axisymmetric_laplacian, bim2a_axisymmetric_advection_diffusion} 
## @end deftypefn

function A = bim1a_axisymmetric_advection_diffusion (x,alpha,gamma,eta,beta)

  ## Check input
  if nargin != 5
    error("bim1a_axisymmetric_advection_diffusion: wrong number of input parameters.");
  elseif !isvector(x)
    error("bim1a_axisymmetric_advection_diffusion: first argument is not a valid vector.");
  endif

  nnodes = length(x);
  nelem  = nnodes-1;

  ## Turn scalar input to a vector of appropriate size
  if isscalar(alpha)
    alpha = alpha*ones(nelem,1);
  endif
  if isscalar(gamma)
    gamma = gamma*ones(nnodes,1);
  endif
  if isscalar(eta)
    eta = eta*ones(nnodes,1);
  endif
  
  if !( isvector(alpha) && isvector(gamma) && isvector(eta) )
    error("bim1a_axisymmetric_advection_diffusion: coefficients are not valid vectors.");
  elseif (length(alpha) != nelem)
    error("bim1a_axisymmetric_advection_diffusion: length of alpha is not equal to the number of elements.");
  elseif (length(gamma) != nnodes)
    error("bim1a_axisymmetric_advection_diffusion: length of gamma is not equal to the number of nodes.");
  elseif (length(eta) != nnodes)
    error("bim1a_axisymmetric_advection_diffusion: length of eta is not equal to the number of nodes.");
  endif

  areak = reshape(diff(x),[],1);
  cm    = reshape((x(1:end-1)+x(2:end))/2,[],1);
  
  if (length(beta) == 1)
    vk = 0;
  elseif (length(beta) == nelem)
    vk = beta .* areak;
  elseif (length(beta) == nnodes)
    vk = diff(beta);
  else
    error("bim1a_axisymmetric_advection_diffusion: coefficient beta has wrong dimensions.");
  endif
  
  gammaetak = bimu_logm ( (gamma.*eta)(1:end-1), (gamma.*eta)(2:end));
  veta      = diff(eta);
  etak      = bimu_logm ( eta(1:end-1), eta(2:end));
  ck        = alpha .* gammaetak .* etak ./ areak .* abs(cm); 

  [bpk, bmk]  = bimu_bernoulli( (vk - veta)./etak);
 
  dm1 = [-(ck.*bmk); NaN]; 
  dp1 = [NaN; -(ck.*bpk)]; 
  d0  = [(ck(1).*bmk(1)); ((ck.*bmk)(2:end) + (ck.*bpk)(1:end-1)); (ck(end).*bpk(end))];
  A   = spdiags([dm1, d0, dp1],-1:1,nnodes,nnodes);

endfunction

%!test
%! n = 3;
%! mesh     = linspace(1,2,n+1)';
%! uex      = @(r) exp(r);
%! duexdr   = @(r) uex(r);
%! d2uexdr2 = @(r) uex(r);
%! Nnodes    = numel(mesh);
%! Nelements = Nnodes-1;
%! D = 1; v = 1;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 1/D*v*ones(Nelements,1);
%! f = @(r) -D./r.*duexdr(r) - D.*d2uexdr2(r) ...
%!         + v./r .* uex(r) + v * duexdr(r);
%! rhs    = bim1a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh));
%! S = bim1a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! u = zeros(Nnodes,1); u([1,end]) = uex(mesh([1 end]));
%! u(2:end-1) = S(2:end-1,2:end-1)\(rhs(2:end-1) - S(2:end-1,[1 end])*u([1 end]));
%! assert(u,uex(mesh),1e-7)

%!test
%! n = 100;
%! mesh      = linspace(0,1,n+1)';
%! cm        = (mesh(1:end-1) + mesh(2:end))/2; 
%! uex       = @(r) - r.^2 + 1;
%! Nnodes    = numel(mesh);
%! Nelements = Nnodes-1;
%! D = 1;  v = 0;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = v;
%! f      = @(r) 4*D;
%! rhs    = bim1a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh));
%! S = bim1a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
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
%! D = 1;  v = cm;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 1/D*v;
%! f      = @(r) 4*D + 2 - 4*r.^2;
%! rhs    = bim1a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh));
%! S = bim1a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! u = zeros(Nnodes,1); u(end) = uex(mesh(end));
%! u(1:end-1) = S(1:end-1,1:end-1)\(rhs(1:end-1) - S(1:end-1,end)*u(end));
%! assert(u,uex(mesh),1e-3)

%!test
%! x = linspace(0,1,101);
%! A = bim1a_axisymmetric_advection_diffusion(x,1,1,1,0);
%! alpha = ones(100,1);
%! gamma = ones(101,1);
%! eta   = gamma;
%! B = bim1a_axisymmetric_advection_diffusion(x,alpha,gamma,eta,0);
%! assert(A,B)
