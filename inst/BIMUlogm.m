function [T] = BIMUlogm(t1,t2)

  ## -*- texinfo -*-
  ## @deftypefn {Function File} @
  ## {[@var{T}]} = BIMUlogm (@var{t1},@var{t2})
  ##
  ## 
  ## 
  ## Input:
  ## @itemize @minus
  ## @item @var{t1}:
  ## @item @var{t2}:
  ## @end itemize
  ##
  ## Output:
  ## @itemize @minus
  ## @item @var{T}:
  ## @end itemize
  ##
  ## @seealso{BIMUbern}
  ## @end deftypefn

  ## This file is part of 
  ##
  ##                   BIM - Box Integration Method Package for Octave
  ##      -------------------------------------------------------------------
  ##              Copyright (C) 2007  Carlo de Falco
  ##              Copyright (C) 2007  Culpo Massimiliano
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
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany
  ##
  ##   Culpo Massimiliano
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany

  T = zeros(size(t2));
  
  sing     = abs(t2-t1)< 100*eps ;
  T(sing)  = (t2(sing)+t1(sing))/2;
  T(~sing) = (t2(~sing)-t1(~sing))./log(t2(~sing)./t1(~sing));
