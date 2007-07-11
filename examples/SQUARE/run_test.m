pkg load msh
pkg load fpl
pkg load bim

## MESH GENERATION AND PROPERTIES
[mesh] = MSH2Mgmsh('quadrato',1,1.2);
[mesh] = BIM2Cmeshproperties(mesh);

## INITIAL GUESS (IF NEEDED)
[xu, yu] = BIM2Cunknowncoord(mesh);
[nelem]  = columns(mesh.t);
[nnodes] = columns(mesh.p);
[uin]    = zeros(nnodes,1);

## COEFFICIENTS FOR THE PROBLEM:
## -div ( \alpha \gamma ( \eta \nabla u - \beta u ) )+ \delta \zeta u = f g
epsilon = .01;
alfa    = ones(nelem,1);
gamma   = ones(nnodes,1);
eta     = epsilon*ones(nnodes,1);
potbeta = xu/epsilon;
delta   = zeros(nelem,1);
zeta    = zeros(nnodes,1);
f       = ones(nnodes,1);
g       = ones(nelem,1);

## ADVECTION AND REACTION MATRIX, RHS
AdvDiff = BIM2Aadvdiff(mesh,alfa,gamma,eta,potbeta);
Mass    = BIM2Areaction(mesh,delta,zeta);
b       = BIM2Arhs(mesh,f,g);

A = AdvDiff + Mass;


## BOUNDARY CONDITIONS
Dlist = BIM2Cunknownsonside(mesh, [2 4]); 	   ## LIST OF DIRICHLET NODES 
Nlist = BIM2Cunknownsonside(mesh, [1 3]); 	   ## LIST OF NEUMANN NODES
Nlist = setdiff(Nlist,Dlist);
Fn    = zeros(length(Nlist),1);           	   ## PRESCRIBED NEUMANN FLUXES
Ilist = setdiff(1:length(uin),union(Dlist,Nlist)); ## LIST OF INTERNAL NODES

## LHS AND RHS PARTITION

## LHS
## DIRICHLET ROW
Add = A(Dlist,Dlist);
Adn = A(Dlist,Nlist);
Adi = A(Dlist,Ilist); 
## NEUMANN ROW
And = A(Nlist,Dlist);
Ann = A(Nlist,Nlist);
Ani = A(Nlist,Ilist); 
## INTERNAL ROW
Aid = A(Ilist,Dlist);
Ain = A(Ilist,Nlist); 
Aii = A(Ilist,Ilist); 

## RHS
bd = b(Dlist);
bn = b(Nlist); 
bi = b(Ilist); 

## UNKNOWNS
ud = uin(Dlist);
un = uin(Nlist); 
ui = uin(Ilist); 

## DISPLACEMENT SOLUTION
temp = [Ann Ani ; Ain Aii ] \ [ Fn+bn-And*ud ; bi-Aid*ud];
un   = temp(1:length(un));
ui   = temp(length(un)+1:end);

## FLUXES THROUGH DIRICHLET SIDES
Fd = Add * ud + Adi * ui + Adn*un - bd;
u  = zeros(nnodes,1);
u(Dlist) = ud; u(Nlist) = un; u(Ilist) = ui;

## GRADIENT OF THE SOLUTION
[gx, gy] = BIM2Cpdegrad(mesh,u);


## FLUXES RECONSTRUCTION
[jxglob,jyglob] = BIM2Cglobalflux(mesh,u,alfa,gamma,eta,potbeta);

## RESULTS VISUALIZATION
FPL2pdesurf(mesh,u);
FPL2pdesurf(mesh,u,"plot_field","gradient");
FPL2quiver(mesh,jxglob,jyglob,"sample_density","1000");
FPL2quiver(mesh,gx,gy,"sample_density","1000");
