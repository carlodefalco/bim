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
##
## @deftypefn {Function File} @
## {@var{A}} = BIM2Aadvdiff (@var{mesh}, @var{alpha}, @var{gamma}, @var{eta}, @var{beta})
##
## Builds the Scharfetter-Gummel matrix for the 
## discretization of the LHS
## of the equation:
## 
## @iftex 
## @tex
## $ -div ( \alpha  \gamma  ( \eta \vect{\nabla} u - \vect{beta} u )) = f $
## @end tex 
## @end iftex 
## @ifinfo
## -div (@var{alpha} * @var{gamma} (@var{eta} grad u - @var{beta} u )) = f
## @end ifinfo
## 
## where: 
## @itemize @minus
## @item @var{alpha}: element-wise constant scalar function
## @item @var{eta}, @var{gamma}: piecewise linear conforming 
## scalar functions
## @item @var{beta}: element-wise constant vector function
## @end itemize
##
## Instead of passing the vector field @var{beta} directly
## one can pass a piecewise linear conforming scalar function
## @var{phi} as the last input.  In such case @var{beta} = grad @var{phi}
## is assumed.  If @var{phi} is a single scalar value @var{beta}
## is assumed to be 0 in the whole domain.
## 
## Example:
## @example
## mesh = MSH2Mstructmesh([0:1/3:1],[0:1/3:1],1,1:4);
## mesh = BIM2Cmeshproperties(mesh);
## x = mesh.p(1,:)';
## Dnodes = BIM2Cunknownsonside(mesh,[2,4]);
## Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
## Varnodes = setdiff(1:Nnodes,Dnodes);
## alpha  = ones(Nelements,1); eta = .1*ones(Nnodes,1);
## beta   = [ones(1,Nelements);zeros(1,Nelements)];
## gamma  = ones(Nnodes,1);
## f      = BIM2Arhs(mesh,ones(Nnodes,1),ones(Nelements,1));
## S = BIM2Aadvdiff(mesh,alpha,gamma,eta,beta);
## u = zeros(Nnodes,1);
## u(Varnodes) = S(Varnodes,Varnodes)\f(Varnodes);
## uex = x - (exp(10*x)-1)/(exp(10)-1);
## assert(u,uex,1e-7)
## @end example
##
## @seealso{BIM2Arhs, BIM2Areaction, BIM2Cmeshproperties}
## @end deftypefn

function [A] = BIM2Aadvdiff(mesh,alpha,gamma,eta,beta)

  Nnodes     = columns(mesh.p);
  Nelements  = columns(mesh.t);
  
  alphaareak = reshape (alpha.*mesh.area,1,1,Nelements);
  shg        = mesh.shg(:,:,:);
  
  
  ## Build local Laplacian matrix	
  
  Lloc = zeros(3,3,Nelements);	
  
  for inode = 1:3
    for jnode = 1:3
      ginode(inode,jnode,:) = mesh.t(inode,:);
      gjnode(inode,jnode,:) = mesh.t(jnode,:);
      Lloc(inode,jnode,:)   = \
	  sum( shg(:,inode,:) .* shg(:,jnode,:),1) .* alphaareak;
    endfor
  endfor
  
  x = mesh.p(1,:);
  x = x(mesh.t(1:3,:));
  y = mesh.p(2,:);
  y = y(mesh.t(1:3,:));
  
  if all(size(beta)==1)
    v12 = 0;
    v23 = 0;
    v31 = 0; 
  elseif all(size(beta)==[2,Nelements])
    v12 = beta(1,:) .* (x(2,:)-x(1,:)) + beta(2,:) .* (y(2,:)-y(1,:));
    v23 = beta(1,:) .* (x(3,:)-x(2,:)) + beta(2,:) .* (y(3,:)-y(2,:));
    v31 = beta(1,:) .* (x(1,:)-x(3,:)) + beta(2,:) .* (y(1,:)-y(3,:)); 
  elseif all(size(beta)==[Nnodes,1])
    betaloc = beta(mesh.t(1:3,:));
    v12     = betaloc(2,:)-betaloc(1,:);
    v23     = betaloc(3,:)-betaloc(2,:);
    v31     = betaloc(1,:)-betaloc(3,:); 
  else
    error("coefficient beta has wrong dimensions");
  endif
  
  etaloc = eta(mesh.t(1:3,:));
  
  eta12 = etaloc(2,:) - etaloc(1,:);
  eta23 = etaloc(3,:) - etaloc(2,:);
  eta31 = etaloc(1,:) - etaloc(3,:);
  
  etalocm1 = BIMUlogm(etaloc(2,:),etaloc(3,:));
  etalocm2 = BIMUlogm(etaloc(3,:),etaloc(1,:));
  etalocm3 = BIMUlogm(etaloc(1,:),etaloc(2,:));
  
  gammaloc = gamma(mesh.t(1:3,:));
  geloc    = gammaloc.*etaloc;
  
  gelocm1 = BIMUlogm(geloc(2,:),geloc(3,:));
  gelocm2 = BIMUlogm(geloc(3,:),geloc(1,:));
  gelocm3 = BIMUlogm(geloc(1,:),geloc(2,:));
  
  [bp12,bm12] = BIMUbern( (v12 - eta12)./etalocm3);
  [bp23,bm23] = BIMUbern( (v23 - eta23)./etalocm1);
  [bp31,bm31] = BIMUbern( (v31 - eta31)./etalocm2);
  
  bp12 = reshape(gelocm3.*etalocm3.*bp12,1,1,Nelements).*Lloc(1,2,:);
  bm12 = reshape(gelocm3.*etalocm3.*bm12,1,1,Nelements).*Lloc(1,2,:);
  bp23 = reshape(gelocm1.*etalocm1.*bp23,1,1,Nelements).*Lloc(2,3,:);
  bm23 = reshape(gelocm1.*etalocm1.*bm23,1,1,Nelements).*Lloc(2,3,:);
  bp31 = reshape(gelocm2.*etalocm2.*bp31,1,1,Nelements).*Lloc(3,1,:);
  bm31 = reshape(gelocm2.*etalocm2.*bm31,1,1,Nelements).*Lloc(3,1,:);
  
  Sloc(1,1,:) = (-bm12-bp31)./reshape(etaloc(1,:),1,1,Nelements);
  Sloc(1,2,:) = bp12./reshape(etaloc(2,:),1,1,Nelements);
  Sloc(1,3,:) = bm31./reshape(etaloc(3,:),1,1,Nelements);
  
  Sloc(2,1,:) = bm12./reshape(etaloc(1,:),1,1,Nelements);
  Sloc(2,2,:) = (-bp12-bm23)./reshape(etaloc(2,:),1,1,Nelements); 
  Sloc(2,3,:) = bp23./reshape(etaloc(3,:),1,1,Nelements);
  
  Sloc(3,1,:) = bp31./reshape(etaloc(1,:),1,1,Nelements);
  Sloc(3,2,:) = bm23./reshape(etaloc(2,:),1,1,Nelements);
  Sloc(3,3,:) = (-bm31-bp23)./reshape(etaloc(3,:),1,1,Nelements);
  
  A = sparse(ginode(:),gjnode(:),Sloc(:));


endfunction

%!test 
%! [mesh] = MSH2Mstructmesh([0:1/3:1],[0:1/3:1],1,1:4);
%! mesh = BIM2Cmeshproperties(mesh);
%! x = mesh.p(1,:)';
%! Dnodes = BIM2Cunknownsonside(mesh,[2,4]);
%! Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
%! Varnodes = setdiff(1:Nnodes,Dnodes);
%! alpha  = ones(Nelements,1); eta = .1*ones(Nnodes,1);
%! beta   = [ones(1,Nelements);zeros(1,Nelements)];
%! gamma  = ones(Nnodes,1);
%! f      = BIM2Arhs(mesh,ones(Nnodes,1),ones(Nelements,1));
%! S = BIM2Aadvdiff(mesh,alpha,gamma,eta,beta);
%! u = zeros(Nnodes,1);
%! u(Varnodes) = S(Varnodes,Varnodes)\f(Varnodes);
%! uex = x - (exp(10*x)-1)/(exp(10)-1);
%! assert(u,uex,1e-7)

%!test 
%! [mesh] = MSH2Mstructmesh([0:1/3:1],[0:1/3:1],1,1:4);
%! mesh = BIM2Cmeshproperties(mesh);
%! x = mesh.p(1,:)';
%! Dnodes = BIM2Cunknownsonside(mesh,[2,4]);
%! Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
%! Varnodes = setdiff(1:Nnodes,Dnodes);
%! alpha  = ones(Nelements,1); eta = .1*ones(Nnodes,1);
%! beta   = x;
%! gamma  = ones(Nnodes,1);
%! f      = BIM2Arhs(mesh,ones(Nnodes,1),ones(Nelements,1));
%! S = BIM2Aadvdiff(mesh,alpha,gamma,eta,beta);
%! u = zeros(Nnodes,1);
%! u(Varnodes) = S(Varnodes,Varnodes)\f(Varnodes);
%! uex = x - (exp(10*x)-1)/(exp(10)-1);
%! assert(u,uex,1e-7)

%!test
%! [mesh] = MSH2Mstructmesh([0:1/3:1],[0:1/3:1],1,1:4);
%! mesh = BIM2Cmeshproperties(mesh);
%! x = mesh.p(1,:)';
%! Dnodes = BIM2Cunknownsonside(mesh,[2,4]);
%! Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
%! Varnodes = setdiff(1:Nnodes,Dnodes);
%! alpha  = 10*ones(Nelements,1); eta = .01*ones(Nnodes,1);
%! beta   = x/10;
%! gamma  = ones(Nnodes,1);
%! f      = BIM2Arhs(mesh,ones(Nnodes,1),ones(Nelements,1));
%! S = BIM2Aadvdiff(mesh,alpha,gamma,eta,beta);
%! u = zeros(Nnodes,1);
%! u(Varnodes) = S(Varnodes,Varnodes)\f(Varnodes);
%! uex = x - (exp(10*x)-1)/(exp(10)-1);
%! assert(u,uex,1e-7)

%!test
%! [mesh] = MSH2Mstructmesh([0:1/3:1],[0:1/3:1],1,1:4);
%! mesh = BIM2Cmeshproperties(mesh);
%! x = mesh.p(1,:)';
%! Dnodes = BIM2Cunknownsonside(mesh,[2,4]);
%! Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
%! Varnodes = setdiff(1:Nnodes,Dnodes);
%! alpha  = 10*ones(Nelements,1); eta = .001*ones(Nnodes,1);
%! beta   = x/100;
%! gamma  = 10*ones(Nnodes,1);
%! f      = BIM2Arhs(mesh,ones(Nnodes,1),ones(Nelements,1));
%! S = BIM2Aadvdiff(mesh,alpha,gamma,eta,beta);
%! u = zeros(Nnodes,1);
%! u(Varnodes) = S(Varnodes,Varnodes)\f(Varnodes);
%! uex = x - (exp(10*x)-1)/(exp(10)-1);
%! assert(u,uex,1e-7)

%!test
%! [mesh] = MSH2Mstructmesh([0:1/1e3:1],[0:1/2:1],1,1:4);
%! mesh = BIM2Cmeshproperties(mesh);
%! x = mesh.p(1,:)';
%! Dnodes = BIM2Cunknownsonside(mesh,[2,4]);
%! Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
%! Varnodes = setdiff(1:Nnodes,Dnodes);
%! alpha  = 3*ones(Nelements,1); eta = x+1;
%! beta   = [ones(1,Nelements);zeros(1,Nelements)];
%! gamma  = 2*x;
%! ff     = 2*(6*x.^2+6*x) - (6*x+6).*(1-2*x)+6*(x-x.^2);
%! f      = BIM2Arhs(mesh,ff,ones(Nelements,1));
%! S = BIM2Aadvdiff(mesh,alpha,gamma,eta,beta);
%! u = zeros(Nnodes,1);
%! u(Varnodes) = S(Varnodes,Varnodes)\f(Varnodes);
%! uex = x - x.^2;
%! assert(u,uex,5e-3)

%!test
%! [mesh] = MSH2Mstructmesh([0:1/1e3:1],[0:1/2:1],1,1:4);
%! mesh = BIM2Cmeshproperties(mesh);
%! x = mesh.p(1,:)';
%! Dnodes = BIM2Cunknownsonside(mesh,[2,4]);
%! Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
%! Varnodes = setdiff(1:Nnodes,Dnodes);
%! alpha  = ones(Nelements,1); eta = ones(Nnodes,1);
%! beta   = 0;
%! gamma  = x+1;
%! ff     = 4*x+1;
%! f      = BIM2Arhs(mesh,ff,ones(Nelements,1));
%! S = BIM2Aadvdiff(mesh,alpha,gamma,eta,beta);
%! u = zeros(Nnodes,1);
%! u(Varnodes) = S(Varnodes,Varnodes)\f(Varnodes);
%! uex = x - x.^2; 
%! assert(u,uex,1e-7)

%!test
%! [mesh] = MSH2Mstructmesh([0:.1:1],[0:.1:1],1,1:4);
%! mesh = BIM2Cmeshproperties(mesh);
%! x = mesh.p(1,:)';y = mesh.p(2,:)';
%! Dnodes = BIM2Cunknownsonside(mesh,[1:4]);
%! Nnodes = columns(mesh.p); Nelements = columns(mesh.t);
%! Varnodes = setdiff(1:Nnodes,Dnodes);
%! alpha  = ones(Nelements,1); diff = 1e-2; eta=diff*ones(Nnodes,1);
%! beta   =[ones(1,Nelements);ones(1,Nelements)];
%! gamma  = x*0+1;
%! ux  = y.*(1-exp((y-1)/diff)) .* (1-exp((x-1)/diff)-x.*exp((x-1)/diff)/diff);
%! uy  = x.*(1-exp((x-1)/diff)) .* (1-exp((y-1)/diff)-y.*exp((y-1)/diff)/diff);
%! uxx = y.*(1-exp((y-1)/diff)) .* (-2*exp((x-1)/diff)/diff-x.*exp((x-1)/diff)/(diff^2));
%! uyy = x.*(1-exp((x-1)/diff)) .* (-2*exp((y-1)/diff)/diff-y.*exp((y-1)/diff)/(diff^2));
%! ff  = -diff*(uxx+uyy)+ux+uy;
%! f   = BIM2Arhs(mesh,ff,ones(Nelements,1));
%! S = BIM2Aadvdiff(mesh,alpha,gamma,eta,beta);
%! u = zeros(Nnodes,1);
%! u(Varnodes) = S(Varnodes,Varnodes)\f(Varnodes);
%! uex = x.*y.*(1-exp((x-1)/diff)).*(1-exp((y-1)/diff)); 
%! assert(u,uex,1e-7)
