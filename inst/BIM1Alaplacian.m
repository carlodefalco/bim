function [A] = BIM1Alaplacian(mesh,epsilon)

  ## -*- texinfo -*-
  ##
  ## @deftypefn {Function File} @
  ## {@var{A}} = BIM1Alaplacian (@var{mesh}, @var{epsilon})
  ##
  ## Builds the finite-element matrix for the 
  ## discretization of the LHS
  ## of the equation:
  ## 
  ## @iftex 
  ## @tex
  ## $ -( \varepsilon  \gamma  ( u' ))' = f $
  ## @end tex 
  ## @end iftex 
  ## @ifinfo
  ## - (@var{epsilon} ( u' ))' = f
  ## @end ifinfo
  ## 
  ## where: 
  ## @itemize @minus
  ## @item @var{epsilon}: elemental values of an piece-wise constant function
  ## @end itemize
  ##
  ##
  ## @seealso{BIM1Arhs, BIM1Areaction}
  ## @end deftypefn

 ## This file is part of 
  ##
  ##                   BIM - Box Integration Method Package for Octave
  ##      -------------------------------------------------------------------
  ##              Copyright (C) 2007  Carlo de Falco and Culpo Massimiliano
  ## 
  ##   BIM is free software; you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation; either version 2 of the License, or
  ##   (at your option) any later version.
  ## 
  ##   BIM is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ## 
  ##   You should have received a copy of the GNU General Public License
  ##   along with BIM; If not, see <http://www.gnu.org/licenses/>.
  ##
  ##
  ##   MAIN AUTHORS:
  ##   Carlo de Falco
  ##   Dublin City University
  ##   Glasnevin, Dublin 9, Ireland
  ##
  ##   Culpo Massimiliano
  ##   Bergische Universitaet Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe fuer Angewandte MathematD-42119 Wuppertal  Gaussstr. 20 
  ##   D-42119 Wuppertal, Germany
    
  Nnodes    = length(mesh);
  h = mesh(2:end)-mesh(1:end-1);
  
  d0 	= [ epsilon(1)./h(1); 
            (epsilon(1:end-1)./h(1:end-1))+(epsilon(2:end)./h(2:end));
            epsilon(end)./h(end)];
  
  d1	= [1000; -epsilon./h];
  dm1	= [ -epsilon./h;1000];
  
  A	= spdiags([dm1, d0, d1],-1:1,Nnodes,Nnodes);