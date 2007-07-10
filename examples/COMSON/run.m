## Mesh generation and properties
[mesh] = MSH2Mgmsh('comson',1,20);
[mesh] = BIM2Cmeshproperties(mesh);

## Construct initial guess
[xu, yu] = BIM2Cunknowncoord(mesh);
nelem = size(mesh.t,2);
nnodi = size(mesh.p,2);
uin   = 0*xu;

## set the coefficients for the problem:
## -div ( \alpha \gamma ( \eta \nabla u - \beta u ) )+ \delta \zeta u = f g
## FIX ME: CAMBIARE I VALORI CON QUALCOSA CHE DIA DEI DISEGNINI PIU' BELLI!
epsilon = .1;
alfa    = ones(nelem,1);
gamma   = ones(nnodi,1);
eta     = epsilon*ones(nnodi,1);
beta    = xu+yu;
delta   = .1*ones(nelem,1);
zeta    = .1*ones(nnodi,1);
f       = ones(nnodi,1);
g       = ones(nelem,1);

## Construction of the discretized operators
AdvDiff = BIM2Aadvdiff(mesh,alfa,gamma,eta,beta);
Mass = BIM2Areaction(mesh,delta,zeta);
b = BIM2Arhs(mesh,f,g);


A = AdvDiff + Mass;


## BOUNDARY CONDITIONS
Dlist = BIM2Cunknownsonside(mesh, [55 38 56 57 35 54,...
				   29 30 1:17 24 25 31 32,...
				   33 34 18:23 26:28]); 	   
## DIRICHLET LIST
Nlist = [] 	   ## NEUMANN LIST
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
FPL2quiver(mesh,jxglob,jyglob,"sample_density","10000");
FPL2quiver(mesh,gx,gy,"sample_density","10000");