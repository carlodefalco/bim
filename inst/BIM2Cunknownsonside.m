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
  ##   along with BIM; if not, write to the Free Software
  ##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  ##   USA
  ##
  ##
  ##   MAIN AUTHOR:
  ##   Carlo de Falco
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany
  ##
  ##   HELP IN REORDERING AND CLEANING THE CODE:
  ##   Culpo Massimiliano
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany

	[nodelist] = MSH2Mnodesonsides(mesh,sidelist);
