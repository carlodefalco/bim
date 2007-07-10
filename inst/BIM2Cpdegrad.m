function [gx, gy] = BIM2Cpdegrad(mesh,u)

	## -*- texinfo -*-
  ##
  ## @deftypefn {Function File} {[@var{gx},@var{gy}]} = BIM2Cpdegrad(@var{mesh},@var{u})
  ##
  ## Builds the P1 approximation of the gradient of the computed solution.
  ##
  ## Input:
  ## @itemize @minus
  ## @item @var{mesh}: PDEtool-like mesh with required field "p", "e", "t".
  ## @item @var{u}: piecewise linear conforming scalar function.
  ## @end itemize 
  ##
  ## @seealso{BIM2Cpdegrad}
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
  ## 	 AID IN PROGRAMMING AND CLEANING THE CODE:
  ##   Culpo Massimiliano
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany




shgx = reshape(mesh.shg(1,:,:),3,[]);
gx = sum(shgx.*u(mesh.t(1:3,:)),1);
shgy = reshape(mesh.shg(2,:,:),3,[]);
gy = sum(shgy.*u(mesh.t(1:3,:)),1);
