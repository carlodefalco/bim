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
## {[@var{bp}, @var{bn}]} = BIMUbern (@var{x})
##
## Computes the values of the Bernoulli function corresponding to x and -x arguments.
## 
## Input:
## @itemize @minus
## @item @var{x}: argument for the Bernoulli function. Could be a matrix of every size.
## @end itemize
##
## Output:
## @itemize @minus
## @item @var{bp}: Bernoulli function for x argument. Same size as @var{x}.
## @item @var{bn}: Bernoulli function for -x argument. Same size as @var{x}.
## @end itemize
##
## @seealso{BIMUlogm}
## @end deftypefn

function [bp,bn] = BIMUbern(x)
  
  xlim= 1e-2;
  ax  = abs(x);
  bp  = zeros(size(x));
  bn  = bp;
  
  block1  = find(~ax);
  block21 = find((ax>80)&x>0);
  block22 = find((ax>80)&x<0);
  block3  = find((ax<=80)&(ax>xlim));
  block4  = find((ax<=xlim)&(ax~=0));
  
  ##  X=0
  bp(block1)=1.;
  bn(block1)=1.;
  
  ## ASYMPTOTICS
  bp(block21)=0.;
  bn(block21)=x(block21);
  bp(block22)=-x(block22);
  bn(block22)=0.;
  
  ## INTERMEDIATE VALUES
  bp(block3)=x(block3)./(exp(x(block3))-1);
  bn(block3)=x(block3)+bp(block3);
  
  ## SMALL VALUES
  if(any(block4))jj=1;
    fp=1.*ones(size(block4));
    fn=fp;
    df=fp;
    segno=1.;
    while (norm(df,inf) > eps),
      jj=jj+1;
      segno=-segno;
      df=df.*x(block4)/jj;
      fp=fp+df;
      fn=fn+segno*df;
    endwhile;
    bp(block4)=1./fp;
    bn(block4)=1./fn;
  endif
  
endfunction
