pkg load msh
pkg load fpl
pkg load bim

## Mesh generation and properties
[mesh] = MSH2Mgmsh('fiume',1,.2);
[mesh] = BIM2Cmeshproperties(mesh);

## Construct initial guess
[xu, yu] = BIM2Cunknowncoord(mesh);
nelem = size(mesh.t,2);
nnodi = size(mesh.p,2);
uin   = 3*xu;

## set the coefficients for the problem:
## -div ( \alpha \gamma ( \eta \nabla u - \beta u ) )+ \delta \zeta u = f g
## FIX ME: CAMBIARE I VALORI CON QUALCOSA CHE DIA DEI DISEGNINI PIU' BELLI!
epsilon = .1;
alfa    = ones(nelem,1);
gamma   = ones(nnodi,1);
eta     = epsilon*ones(nnodi,1);
beta    = xu+yu;
delta   = ones(nelem,1);
zeta    = ones(nnodi,1);
f       = ones(nnodi,1);
g       = ones(nelem,1);

## Construction of the discretized operators
AdvDiff = BIM2Aadvdiff(mesh,alfa,gamma,eta,beta);
Mass = BIM2Areaction(mesh,delta,zeta);
b = BIM2Arhs(mesh,f,g);


A = AdvDiff + Mass;


## BOUNDARY CONDITIONS
Dlist = BIM2Cunknownsonside(mesh, [8 18]); 	   ## DIRICHLET LIST
Nlist = BIM2Cunknownsonside(mesh, [23 24]); 	   ## NEUMANN LIST
Nlist = setdiff(Nlist,Dlist);
Fn    = zeros(length(Nlist),1);           	   ## PRESCRIBED NEUMANN FLUXES
Ilist = setdiff(1:length(uin),union(Dlist,Nlist)); ## INTERNAL LIST

## Partition the LHS and RHS

Add = A(Dlist,Dlist);
Adn = A(Dlist,Nlist); ## shoud be all zeros hopefully!!
Adi = A(Dlist,Ilist); 

And = A(Nlist,Dlist); ## shoud be all zeros hopefully!!
Ann = A(Nlist,Nlist);
Ani = A(Nlist,Ilist); 

Aid = A(Ilist,Dlist);
Ain = A(Ilist,Nlist); 
Aii = A(Ilist,Ilist); 

bd = b(Dlist);
bn = b(Nlist); 
bi = b(Ilist); 

ud = uin(Dlist);
un = uin(Nlist); 
ui = uin(Ilist); 

## SOLUTION
temp = [Ann Ani ; Ain Aii ] \ [ Fn+bn-And*ud ; bi-Aid*ud];
un   = temp(1:length(un));
ui   = temp(length(un)+1:end);

## Compute fluxes through Dirichlet sides

Fd = Add * ud + Adi * ui + Adn*un - bd;

u = zeros(nnodi,1);
u(Dlist) = ud; u(Nlist) = un; u(Ilist) = ui;

## Reconstruction of the currents

## GRADIENT OF THE SOLUTION
[gx, gy] = BIM2Cpdegrad(mesh,u);


## FLUXES
[jxglob,jyglob] = BIM2Cglobalflux(mesh,u,alfa,gamma,eta,beta);

## VISUALIZATION
FPL2pdesurf(mesh,u);
FPL2quiver(mesh,jxglob,jyglob,"sample_density","1000");
FPL2quiver(mesh,gx,gy,"sample_density","1000");
