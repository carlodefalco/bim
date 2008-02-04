function b = BIM2Arhs(mesh,f,g)

  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{b}]} = BIM2Arhs(@var{mesh}, @var{f}, @var{g})
  ##
  ## Constructs the RHS for the DAR problem:
  ## @iftex 
  ## @tex
  ## $ -div ( \alpha  \gamma  ( \eta \vect{\nabla} u - \vect{beta} u )) + \delta \zeta u = f g $
  ## @end tex 
  ## @end iftex 
  ## @ifinfo
  ## -div (@var{alpha} * @var{gamma} (@var{eta} grad u - @var{beta} u )) + @var{delta} * @var{zeta} u = @var{f}*@var{g}
  ## @end ifinfo
  ## 
  ## Input:
  ## @itemize @minus
  ## @item @var{mesh}: PDEtool-like mesh with required field "p", "e", "t".
  ## @item @var{g}: element-wise constant scalar function.
  ## @item @var{f}: piecewise linear conforming scalar function.
  ## @end itemize 
  ##
  ## @seealso{BIM2Areaction, BIM2Aadvdiff, BIM2Cmeshproperties}
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
  
  Nnodes    =size(mesh.p,2);
  Nelements =size(mesh.t,2);
  
  f       = f(mesh.t(1:3,:));
  wjacdet = mesh.wjacdet;
  %% build local matrix	
  Blocmat=zeros(3,Nelements);	
  for inode=1:3
    Blocmat(inode,:) = g'.*f(inode,:).*wjacdet(inode,:);
  end

  gnode=(mesh.t(1:3,:));
  %% assemble global matrix

  b = sparse(gnode(:),1,Blocmat(:));
