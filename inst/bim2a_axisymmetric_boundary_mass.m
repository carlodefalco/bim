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
## @deftypefn {Function File} {[@var{M}]} = @
## bim2a_axisymmetric_boundary_mass(@var{mesh},@var{sidelist},@var{nodelist})
##
## Build the lumped boundary mass matrix needed to apply Robin and Neumann
## boundary conditions in a problem in cylindrical coordinates with 
## axisymmetric configuration.
##
## The vector @var{sidelist} contains the list of the side edges
## contributing to the mass matrix.
##
## The optional argument @var{nodelist} contains the list of the
## degrees of freedom on the boundary.
##
## @seealso{bim2a_axisymmetric_rhs, bim2a_axisymmetric_advection_diffusion, 
## bim2a_axisymmetric_laplacian, bim2a_axisymmetric_reaction, bim2a_boundary_mass} 
## @end deftypefn

function [M] = bim2a_axisymmetric_boundary_mass(mesh,sidelist,nodelist)

  ## Check input
  if (nargin > 3)
    error ("bim2a_axisymmetric_boundary_mass: wrong number of input parameters.");
  elseif !(isstruct(mesh)     && isfield(mesh,"p") &&
	       isfield (mesh,"t") && isfield(mesh,"e"))
    error("bim2a_axisymmetric_boundary_mass: first input is not a valid mesh structure.");
  elseif !( isvector(sidelist) && isnumeric(sidelist) )
    error("bim2a_axisymmetric_boundary_mass: second input is not a valid numeric vector.");
  elseif !(all(mesh.p(1,:) >= 0) || all(mesh.p(1,:) <= 0))
    error("bim2a_axisymmetric_boundary_mass: the input mesh cannot intersect the rotation axis r=0.");
  endif

  if (nargin < 3)
    [nodelist] = bim2c_unknowns_on_side(mesh,sidelist);
  endif

  r = abs (mesh.p(1,nodelist));

  edges = [];
  for ie = sidelist
    edges = [ edges, mesh.e([1:2 5],mesh.e(5,:)==ie)];
  endfor
  l  = sqrt((mesh.p(1,edges(1,:))-mesh.p(1,edges(2,:))).^2 +
	          (mesh.p(2,edges(1,:))-mesh.p(2,edges(2,:))).^2);
  
  dd = zeros(size(nodelist));
  
  for in = 1:length(nodelist)
    dd (in) = ( sum(r(in).*l(edges(1,:)==nodelist(in))) ...
              + sum(r(in).*l(edges(2,:)==nodelist(in))) )/2;
  endfor
  
  M = sparse(diag(dd));

endfunction

%!test
%! n = 3;
%! [mesh] = msh2m_structured_mesh(linspace(1,2,n+1),linspace(0,1,n+1),1,1:4);
%! mesh   = bim2c_mesh_properties(mesh);
%! uex    = @(r,z) exp(r);
%! duexdr = @(r,z) uex(r,z);
%! d2uexdr2 = @(r,z) uex(r,z);
%! duexdz = @(r,z) 0*uex(r,z);
%! d2uexdz2 = @(r,z) 0*uex(r,z);
%! Rnodesr  = bim2c_unknowns_on_side(mesh,[2]);
%! Rnodesl  = bim2c_unknowns_on_side(mesh,[4]);
%! Nnodes    = columns(mesh.p);
%! Nelements = columns(mesh.t);
%! D = 1; vr = 1; vz = 0;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 1/D*[vr*ones(1,Nelements); vz*ones(1,Nelements)];
%! f = @(r,z) -D./r.*duexdr(r,z) - D.*d2uexdr2(r,z) ...
%!           + vr./r .* uex(r,z) + vr * duexdr(r,z) ...
%!           - D.*d2uexdz2(r,z) + vz * duexdz(r,z);
%! gr = @(r,z) uex(r,z) - 1 * (-D*duexdr(r,z) + vr*uex(r,z));
%! gl = @(r,z) uex(r,z) - (-1) * (-D*duexdr(r,z) + vr*uex(r,z));
%! rhs    = bim2a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh.p(1,:), mesh.p(2,:)));
%! S = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! Mr = bim2a_axisymmetric_boundary_mass(mesh,2);  Ml = bim2a_axisymmetric_boundary_mass(mesh,4);
%! S(Rnodesr,Rnodesr) += Mr;
%! rhs(Rnodesr) += diag(Mr) .* gr(mesh.p(1,Rnodesr), mesh.p(2,Rnodesr))';
%! S(Rnodesl,Rnodesl) += Ml;
%! rhs(Rnodesl) += diag(Ml) .* gl(mesh.p(1,Rnodesl), mesh.p(2,Rnodesl))';
%! u = S\rhs; 
%! assert(u,uex(mesh.p(1,:), mesh.p(2,:))',1e-7)

%!test
%! n = 10;
%! [mesh] = msh2m_structured_mesh(linspace(1,2,n+1),linspace(0,1,n+1),1,1:4);
%! mesh   = bim2c_mesh_properties(mesh);
%! uex    = @(r,z) exp(r);
%! duexdr = @(r,z) uex(r,z);
%! d2uexdr2 = @(r,z) uex(r,z);
%! duexdz = @(r,z) 0*uex(r,z);
%! d2uexdz2 = @(r,z) 0*uex(r,z);
%! Rnodesr  = bim2c_unknowns_on_side(mesh,[2]);
%! Rnodesl  = bim2c_unknowns_on_side(mesh,[4]);
%! Rnodesb  = bim2c_unknowns_on_side(mesh,[1]);
%! Rnodest  = bim2c_unknowns_on_side(mesh,[3]);
%! Nnodes    = columns(mesh.p);
%! Nelements = columns(mesh.t);
%! D = 1; vr = 1; vz = 0;
%! alpha  = D*ones(Nelements,1);
%! gamma  = ones(Nnodes,1);
%! eta    = ones(Nnodes,1);
%! beta   = 1/D*[vr*ones(1,Nelements); vz*ones(1,Nelements)];
%! f = @(r,z) -D./r.*duexdr(r,z) - D.*d2uexdr2(r,z) ...
%!           + vr./r .* uex(r,z) + vr * duexdr(r,z) ...
%!           - D.*d2uexdz2(r,z) + vz * duexdz(r,z);
%! gr = @(r,z) uex(r,z) - 1 * (-D*duexdr(r,z) + vr*uex(r,z));
%! gl = @(r,z) uex(r,z) - (-1) * (-D*duexdr(r,z) + vr*uex(r,z));
%! gb = @(r,z) uex(r,z) - (-1) * (-D*duexdz(r,z) + vz*uex(r,z));
%! gt = @(r,z) uex(r,z) - 1 * (-D*duexdz(r,z) + vz*uex(r,z));
%! rhs    = bim2a_axisymmetric_rhs(mesh, ones(Nelements,1), f(mesh.p(1,:), mesh.p(2,:)));
%! S = bim2a_axisymmetric_advection_diffusion(mesh,alpha,gamma,eta,beta);
%! Mr = bim2a_axisymmetric_boundary_mass(mesh,2);  Ml = bim2a_axisymmetric_boundary_mass(mesh,4);
%! Mb = bim2a_axisymmetric_boundary_mass(mesh,1);  Mt = bim2a_axisymmetric_boundary_mass(mesh,3);
%! S(Rnodesr,Rnodesr) += Mr;
%! rhs(Rnodesr) += diag(Mr) .* gr(mesh.p(1,Rnodesr), mesh.p(2,Rnodesr))';
%! S(Rnodesl,Rnodesl) += Ml;
%! rhs(Rnodesl) += diag(Ml) .* gl(mesh.p(1,Rnodesl), mesh.p(2,Rnodesl))';
%! S(Rnodesb,Rnodesb) += Mb;
%! rhs(Rnodesb) += diag(Mb) .* gb(mesh.p(1,Rnodesb), mesh.p(2,Rnodesb))';
%! S(Rnodest,Rnodest) += Mt;
%! rhs(Rnodest) += diag(Mt) .* gt(mesh.p(1,Rnodest), mesh.p(2,Rnodest))';
%! u = S\rhs; 
%! assert(u,uex(mesh.p(1,:), mesh.p(2,:))',1e-7)
