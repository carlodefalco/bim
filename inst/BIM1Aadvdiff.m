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
## {@var{A}} = BIM1Aadvdiff (@var{mesh}, @var{alpha}, @var{gamma}, @var{eta}, @var{beta})
##
## Builds the Scharfetter-Gummel matrix for the 
## discretization of the LHS
## of the equation:
## @iftex 
## @tex
## $ - ( \alpha  \gamma  ( \eta  u' - \vect{beta} u ))' = f $
## @end tex 
## @end iftex 
## @ifinfo
## - (@var{alpha} * @var{gamma} (@var{eta} u' - @var{beta} u ))' = f
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
##
## @seealso{BIM1Arhs, BIM21reaction, BIM1Alaplacian}
## @end deftypefn

function A = BIM1Aadvdiff(x,alpha,gamma,eta,beta)

  Nnodes     = length(x);
  Nelements  = Nnodes-1;
  
  areak = diff(x);
 
  if (length(beta) == 1)
    vk = 0;
  elseif (length(beta) == Nelements)
    vk = beta .* areak;
  elseif (length(beta) == Nnodes)
    vk = diff(beta);
  else
    error("coefficient beta has wrong dimensions");
  endif
  
  gammaetak = BIMUlogm ( (gamma.*eta)(1:end-1), (gamma.*eta)(2:end));
  veta      = diff(eta);
  etak      = BIMUlogm ( eta(1:end-1), eta(2:end));
  ck        = alpha .* gammaetak .* etak ./ areak; 

  [bpk, bmk]  = BIMUbern( (vk - veta)./etak);
 
  dm1 = [-(ck.*bmk); NaN]; 
  dp1 = [NaN; -(ck.*bpk)]; 
  d0  = [(ck(1).*bpk(1)); ((ck.*bmk)(2:end) + (ck.*bpk)(1:end-1)); (ck(end).*bmk(end))];
  A   = spdiags([dm1, d0, dp1],-1:1,Nnodes,Nnodes);

endfunction


