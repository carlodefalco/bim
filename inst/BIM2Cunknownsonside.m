function [nodelist] = BIM2Cunknownsonside(mesh, sidelist)
	
  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{nodelist}]} = BIM2Cunknownsonside(@var{mesh},@var{sidelist})
  ##
  ## Returns the list of the mesh nodes that lie on the specified geometrical sides.
  ##
  ## Input:
  ## @itemize @minus
  ## @item @var{mesh}: PDEtool-like mesh with required field "p", "e", "t".
  ## @item @var{sidelist}: list of the sides of the geometrical border.
  ## @end itemize 
  ##
  ## Output:
  ## @itemize @minus
  ## @item @var{nodelist}: list of the nodes that lie on the specified sides.
  ## @end itemize
  ##
  ## @seealso{BIM2Cunknowncoord}
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

  [nodelist] = MSH2Mnodesonsides(mesh,sidelist);
