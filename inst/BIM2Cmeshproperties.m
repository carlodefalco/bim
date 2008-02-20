function [omesh] = BIM2Cmeshproperties(imesh)

  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{omesh}]} = BIM2Cmeshproperties(@var{imesh})
  ##
  ## Creates an omesh structure starting from imesh. All the properties needed by BIM are added as fields.
  ##
  ## Input:
  ## @itemize @minus
  ## @item @var{imesh}: PDEtool-like mesh with required field "p", "e", "t".
  ## @end itemize 
  ##
  ## Output:
  ## @itemize @minus
  ## @item @var{omesh}: PDEtool-like mesh structure with added fields needed by BIM method.
  ## @end itemize
  ##
  ## @seealso{BIM2Areaction, BIM2Aadvdiff, BIM2Crhs}
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
  
  omesh = imesh;
  [omesh.wjacdet,omesh.area,omesh.shg] = MSH2Mgeomprop(imesh,"wjacdet","area","shg");
