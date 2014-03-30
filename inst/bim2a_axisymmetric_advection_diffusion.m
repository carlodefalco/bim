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
## @deftypefn {Function File} @
## {[@var{A}]} = @
## bim2a_axisymmetric_advection_diffusion(@var{mesh},@var{alpha},@var{gamma},@var{eta},@var{beta})
##
## Build the Scharfetter-Gummel stabilized stiffness matrix for a
## diffusion-advection problem in cylindrical coordinates with axisymmetric 
## configuration. Rotational symmetry is assumed with respect to be the vertical
## axis r=0. Only plane geometries that DO NOT intersect the symmetry axis 
## are admitted.
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
## F = [Fr, Fz]' = - @var{alpha} * @var{gamma} ( @var{eta} grad (u) - @var{beta} u )
##
## where @var{alpha} is an element-wise constant scalar function,
## @var{eta} and @var{gamma} are piecewise linear conforming scalar
## functions, @var{beta} is an element-wise constant vector function.
##
## Instead of passing the vector field @var{beta} directly, one can pass
## a piecewise linear conforming scalar function  @var{phi} as the last
## input. In such case @var{beta} = grad @var{phi} is assumed.
##
## If @var{phi} is a single scalar value @var{beta} is assumed to be 0
## in the whole domain.
##
## @seealso{bim2a_axisymmetric_rhs, bim2a_axisymmetric_reaction, 
## bim2a_advection_diffusion, bim2c_mesh_properties}
## @end deftypefn

function [A] = bim2a_axisymmetric_advection_diffusion (mesh, alpha, gamma, eta, beta)

  ## Check input
  if nargin != 5
    error("bim2a_axisymmetric_advection_diffusion: wrong number of input parameters.");
  elseif !(isstruct(mesh) && isfield(mesh,"p") &&
	         isfield (mesh,"t") && isfield(mesh,"e"))
    error("bim2a_axisymmetric_advection_diffusion: first input is not a valid mesh structure.");
  elseif !(all(mesh.p(1,:) >= 0) || all(mesh.p(1,:) <= 0))
    error("bim2a_axisymmetric_advection_diffusion: the input mesh cannot intersect the rotation axis r=0.");
  endif
  
  nnodes = columns(mesh.p);
  nelem  = columns(mesh.t);

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
    error("bim2a_axisymmetric_advection_diffusion: coefficients are not valid vectors.");
  elseif length(alpha) != nelem
    error("bim2a_axisymmetric_advection_diffusion: length of alpha is not equal to the number of elements.");
  elseif length(gamma) != nnodes
    error("bim2a_axisymmetric_advection_diffusion: length of gamma is not equal to the number of nodes.");
  elseif length(eta) != nnodes
    error("bim2a_axisymmetric_advection_diffusion: length of eta is not equal to the number of nodes.");
  endif
  
  x = abs( mesh.p(1,:));
  x = x(mesh.t(1:3,:));
  y = mesh.p(2,:);
  y = y(mesh.t(1:3,:));

  rcm = sum (x, 1) / (rows (mesh.t) - 1);
   
  alphaareak = reshape (alpha.*mesh.area.*rcm',1,1,nelem);
  shg        = mesh.shg(:,:,:);
  
  ## Build local Laplacian matrix
  Lloc = zeros(3,3,nelem);	
  
  for inode = 1:3
    for jnode = 1:3
      ginode(inode,jnode,:) = mesh.t(inode,:);
      gjnode(inode,jnode,:) = mesh.t(jnode,:);
      Lloc(inode,jnode,:)   = sum( shg(:,inode,:) .* shg(:,jnode,:),1) .* alphaareak;
    endfor
  endfor
    
  if all(size(beta)==1)
    v12 = 0;
    v23 = 0;
    v31 = 0; 
  elseif all(size(beta)==[2,nelem])
    v12 = beta(1,:) .* (x(2,:)-x(1,:)) + beta(2,:) .* (y(2,:)-y(1,:));
    v23 = beta(1,:) .* (x(3,:)-x(2,:)) + beta(2,:) .* (y(3,:)-y(2,:));
    v31 = beta(1,:) .* (x(1,:)-x(3,:)) + beta(2,:) .* (y(1,:)-y(3,:)); 
  elseif all(size(beta)==[nnodes,1])
    betaloc = beta(mesh.t(1:3,:));
    v12     = betaloc(2,:)-betaloc(1,:);
    v23     = betaloc(3,:)-betaloc(2,:);
    v31     = betaloc(1,:)-betaloc(3,:); 
  else
    error("bim2a_axisymmetric_advection_diffusion: coefficient beta has wrong dimensions.");
  endif
  
  etaloc = eta(mesh.t(1:3,:));
  
  eta12 = etaloc(2,:) - etaloc(1,:);
  eta23 = etaloc(3,:) - etaloc(2,:);
  eta31 = etaloc(1,:) - etaloc(3,:);
  
  etalocm1 = bimu_logm(etaloc(2,:),etaloc(3,:));
  etalocm2 = bimu_logm(etaloc(3,:),etaloc(1,:));
  etalocm3 = bimu_logm(etaloc(1,:),etaloc(2,:));
  
  gammaloc = gamma(mesh.t(1:3,:));
  geloc    = gammaloc.*etaloc;
  
  gelocm1 = bimu_logm (geloc(2,:), geloc(3,:));
  gelocm2 = bimu_logm (geloc(3,:), geloc(1,:));
  gelocm3 = bimu_logm (geloc(1,:), geloc(2,:));
  
  [bp12,bm12] = bimu_bernoulli ((v12 - eta12) ./ etalocm3);
  [bp23,bm23] = bimu_bernoulli ((v23 - eta23) ./ etalocm1);
  [bp31,bm31] = bimu_bernoulli ((v31 - eta31) ./ etalocm2);
  
  bp12 = reshape(gelocm3.*etalocm3.*bp12,1,1,nelem).*Lloc(1,2,:);
  bm12 = reshape(gelocm3.*etalocm3.*bm12,1,1,nelem).*Lloc(1,2,:);
  bp23 = reshape(gelocm1.*etalocm1.*bp23,1,1,nelem).*Lloc(2,3,:);
  bm23 = reshape(gelocm1.*etalocm1.*bm23,1,1,nelem).*Lloc(2,3,:);
  bp31 = reshape(gelocm2.*etalocm2.*bp31,1,1,nelem).*Lloc(3,1,:);
  bm31 = reshape(gelocm2.*etalocm2.*bm31,1,1,nelem).*Lloc(3,1,:);
  
  Sloc(1,1,:) = (-bm12-bp31)./reshape(etaloc(1,:),1,1,nelem);
  Sloc(1,2,:) = bp12./reshape(etaloc(2,:),1,1,nelem);
  Sloc(1,3,:) = bm31./reshape(etaloc(3,:),1,1,nelem);
  
  Sloc(2,1,:) = bm12./reshape(etaloc(1,:),1,1,nelem);
  Sloc(2,2,:) = (-bp12-bm23)./reshape(etaloc(2,:),1,1,nelem); 
  Sloc(2,3,:) = bp23./reshape(etaloc(3,:),1,1,nelem);
  
  Sloc(3,1,:) = bp31./reshape(etaloc(1,:),1,1,nelem);
  Sloc(3,2,:) = bm23./reshape(etaloc(2,:),1,1,nelem);
  Sloc(3,3,:) = (-bm31-bp23)./reshape(etaloc(3,:),1,1,nelem);
  
  A = sparse(ginode(:),gjnode(:),Sloc(:));

endfunction

%!test
%! n         = 3;
%! [mesh]    = msh2m_structured_mesh(linspace(1,2,n+1),linspace(0,1,n+1),1,1:4);
%! mesh      = bim2c_mesh_properties(mesh);
%! uex       = @(r,z) exp(r);
%! duexdr    = @(r,z) uex(r,z);
%! d2uexdr2  = @(r,z) uex(r,z);
%! duexdz    = @(r,z) 0*uex(r,z);
%! d2uexdz2  = @(r,z) 0*uex(r,z);
%! Dnodes    = bim2c_unknowns_on_side(mesh,[2,4]);
%! Nnodes    = columns(mesh.p);
%! Nelements = columns(mesh.t);
%! Varnodes  = setdiff(1:Nnodes,Dnodes);
%! D = 1; vr = 1; vz = 0;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 1/D*[vr*ones(1,Nelements); vz*ones(1,Nelements)];
%! f = @(r,z) -D./r.*duexdr(r,z) - D.*d2uexdr2(r,z) ...
%!           + vr./r .* uex(r,z) + vr * duexdr(r,z) ...
%!           - D.*d2uexdz2(r,z) + vz * duexdz(r,z);
%! rhs = bim2a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh.p(1,:), mesh.p(2,:)));
%! S   = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! u   = zeros(Nnodes,1); u(Dnodes) = uex(mesh.p(1,Dnodes), mesh.p(2,Dnodes));
%! u(Varnodes) = S(Varnodes,Varnodes)\(rhs(Varnodes) - S(Varnodes,Dnodes)*u(Dnodes));
%! assert(u,uex(mesh.p(1,:), mesh.p(2,:))',1e-7)

%!test
%! n         = 20;
%! [mesh]    = msh2m_structured_mesh(linspace(1,2,n+1),linspace(0,1,n+1),1,1:4);
%! mesh      = bim2c_mesh_properties(mesh);
%! uex       = @(r,z) exp(r) .* exp(1-z);
%! duexdr    = @(r,z) uex(r,z);
%! d2uexdr2  = @(r,z) uex(r,z);
%! duexdz    = @(r,z) -uex(r,z);
%! d2uexdz2  = @(r,z) uex(r,z);
%! Dnodes    = bim2c_unknowns_on_side(mesh,[1,2,3,4]);
%! Nnodes    = columns(mesh.p);
%! Nelements = columns(mesh.t);
%! Varnodes  = setdiff(1:Nnodes,Dnodes);
%! D = 1; vr = 1; vz = 1;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 1/D*[vr*ones(1,Nelements); vz*ones(1,Nelements)];
%! f = @(r,z) -D./r.*duexdr(r,z) - D.*d2uexdr2(r,z) ...
%!           + vr./r .* uex(r,z) + vr * duexdr(r,z) ...
%!           - D.*d2uexdz2(r,z) + vz * duexdz(r,z);
%! rhs = bim2a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh.p(1,:), mesh.p(2,:)));
%! S   = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! u   = zeros(Nnodes,1); u(Dnodes) = uex(mesh.p(1,Dnodes), mesh.p(2,Dnodes));
%! u(Varnodes) = S(Varnodes,Varnodes)\(rhs(Varnodes) - S(Varnodes,Dnodes)*u(Dnodes));
%! assert(u,uex(mesh.p(1,:), mesh.p(2,:))',1e-3)

%!test
%! n         = 10;
%! [mesh]    = msh2m_structured_mesh(linspace(1,2,n+1),linspace(0,1,n+1),1,1:4);
%! mesh      = bim2c_mesh_properties(mesh);
%! uex       = @(r,z) exp(r) .* exp(1-z);
%! duexdr    = @(r,z) uex(r,z);
%! d2uexdr2  = @(r,z) uex(r,z);
%! duexdz    = @(r,z) -uex(r,z);
%! d2uexdz2  = @(r,z) uex(r,z);
%! Dnodes    = bim2c_unknowns_on_side(mesh,[1,2,3,4]);
%! Nnodes    = columns(mesh.p);
%! Nelements = columns(mesh.t);
%! Varnodes  = setdiff(1:Nnodes,Dnodes);
%! D      = 1;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 1/D * mesh.p(1,:)';
%! f = @(r,z) -D./r.*duexdr(r,z) - D.*d2uexdr2(r,z) ...
%!           + 1./r .* uex(r,z) + duexdr(r,z) ...
%!           - D.*d2uexdz2(r,z);
%! rhs = bim2a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh.p(1,:), mesh.p(2,:)));
%! S   = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! u   = zeros(Nnodes,1); u(Dnodes) = uex(mesh.p(1,Dnodes), mesh.p(2,Dnodes));
%! u(Varnodes) = S(Varnodes,Varnodes)\(rhs(Varnodes) - S(Varnodes,Dnodes)*u(Dnodes));
%! assert(u,uex(mesh.p(1,:), mesh.p(2,:))',1e-3)

%!test
%! n         = 10;
%! [mesh]    = msh2m_structured_mesh(linspace(1,2,n+1),linspace(0,1,n+1),1,1:4);
%! mesh      = bim2c_mesh_properties(mesh);
%! uex       = @(r,z) exp(r) .* exp(1-z);
%! duexdr    = @(r,z) uex(r,z);
%! d2uexdr2  = @(r,z) uex(r,z);
%! duexdz    = @(r,z) -uex(r,z);
%! d2uexdz2  = @(r,z) uex(r,z);
%! Dnodes    = bim2c_unknowns_on_side(mesh,[1,2,3,4]);
%! Nnodes    = columns(mesh.p);
%! Nelements = columns(mesh.t);
%! Varnodes  = setdiff(1:Nnodes,Dnodes);
%! D         = 1;
%! alpha     = D*ones(Nelements,1);
%! gamma     = ones(Nnodes,1);
%! eta       = ones(Nnodes,1);
%! beta      = 1/D * 1/2*(mesh.p(1,:)').^2;
%! f         = @(r,z) 1./r.*(1+r).*(r-D) .* uex(r,z);
%! rhs = bim2a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh.p(1,:), mesh.p(2,:)));
%! S   = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! u   = zeros(Nnodes,1); u(Dnodes) = uex(mesh.p(1,Dnodes), mesh.p(2,Dnodes));
%! u(Varnodes) = S(Varnodes,Varnodes)\(rhs(Varnodes) - S(Varnodes,Dnodes)*u(Dnodes));
%! assert(u,uex(mesh.p(1,:), mesh.p(2,:))',1e-3)

%!test
%! n         = 3;
%! [mesh]    = msh2m_structured_mesh(linspace(-2,-1,n+1),linspace(0,1,n+1),1,1:4);
%! mesh      = bim2c_mesh_properties(mesh);
%! uex       = @(r,z) exp(r);
%! duexdr    = @(r,z) uex(r,z);
%! d2uexdr2  = @(r,z) uex(r,z);
%! duexdz    = @(r,z) 0*uex(r,z);
%! d2uexdz2  = @(r,z) 0*uex(r,z);
%! Dnodes    = bim2c_unknowns_on_side(mesh,[2,4]);
%! Nnodes    = columns(mesh.p);
%! Nelements = columns(mesh.t);
%! Varnodes  = setdiff(1:Nnodes,Dnodes);
%! D = 1; vr = 1; vz = 0;
%! alpha     = D*ones(Nelements,1);
%! gamma     = ones(Nnodes,1);
%! eta       = ones(Nnodes,1);
%! beta      = 1/D*[vr*ones(1,Nelements); vz*ones(1,Nelements)];
%! f = @(r,z) -D./r.*duexdr(r,z) - D.*d2uexdr2(r,z) ...
%!           + vr./r .* uex(r,z) + vr * duexdr(r,z) ...
%!           - D.*d2uexdz2(r,z) + vz * duexdz(r,z);
%! rhs = bim2a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh.p(1,:), mesh.p(2,:)));
%! S   = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! u   = zeros(Nnodes,1); u(Dnodes) = uex(mesh.p(1,Dnodes), mesh.p(2,Dnodes));
%! u(Varnodes) = S(Varnodes,Varnodes)\(rhs(Varnodes) - S(Varnodes,Dnodes)*u(Dnodes));
%! assert(u,uex(mesh.p(1,:), mesh.p(2,:))',1e-7)

%!test
%! n         = 10;
%! [mesh]    = msh2m_structured_mesh(linspace(-2,-1,n+1),linspace(0,1,n+1),1,1:4);
%! mesh      = bim2c_mesh_properties(mesh);
%! uex       = @(r,z) exp(r) .* exp(1-z);
%! duexdr    = @(r,z) uex(r,z);
%! d2uexdr2  = @(r,z) uex(r,z);
%! duexdz    = @(r,z) -uex(r,z);
%! d2uexdz2  = @(r,z) uex(r,z);
%! Dnodes    = bim2c_unknowns_on_side(mesh,[1,2,3,4]);
%! Nnodes    = columns(mesh.p);
%! Nelements = columns(mesh.t);
%! Varnodes  = setdiff(1:Nnodes,Dnodes);
%! D         = 1;
%! alpha     = D*ones(Nelements,1);
%! gamma     = ones(Nnodes,1);
%! eta       = ones(Nnodes,1);
%! beta      = 1/D * 1/2*(mesh.p(1,:)').^2;
%! f         = @(r,z) 1./r.*(1+r).*(r-D) .* uex(r,z);
%! rhs = bim2a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh.p(1,:), mesh.p(2,:)));
%! S   = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! u   = zeros(Nnodes,1); u(Dnodes) = uex(mesh.p(1,Dnodes), mesh.p(2,Dnodes));
%! u(Varnodes) = S(Varnodes,Varnodes)\(rhs(Varnodes) - S(Varnodes,Dnodes)*u(Dnodes));
%! assert(u,uex(mesh.p(1,:), mesh.p(2,:))',1e-3)

%!test
%! [mesh] = msh2m_structured_mesh([0:.1:1],[0:.1:1],1,1:4);
%! mesh   = bim2c_mesh_properties(mesh);
%! x = mesh.p(1,:)'; y = mesh.p(2,:)';
%! Dnodes = bim2c_unknowns_on_side(mesh,[1:4]);
%! Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
%! alpha  = ones(Nelements,1); eta=ones(Nnodes,1);
%! beta   = 0;
%! gamma  = ones(Nnodes,1);
%! A = bim2a_axisymmetric_advection_diffusion(mesh,1,1,1,0);
%! B = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! assert(A,B)
