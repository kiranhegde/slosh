module par
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    !boundary condition
    integer::NR,NZ,Nout,number_write
    integer::RS,RE,ZS,ZE,BCXLEF,BCXRIG,BCYBOT,BCYTOP
    integer::set_compute

    integer::NXD,NZD
    real(kind=dp)::residual,residualold
    real(kind=dp)::mass   !1.e-8_dp  !residual
    real(kind=dp)::thetaA,thetaM
    real(kind=dp)::epiron 
    real(kind=dp)::h0,r0           !height
    real(kind=dp)::deg          !angle of free surface
    real(kind=dp)::XA,XB,YA,YB,dr,dz
    
    real(kind=dp),parameter::pi=4_dp*atan(1.)
    integer::i,k,fk
    real(kind=dp),pointer,dimension(:)::RP,RC,ZP,ZC
    !transform
    real(kind=dp),pointer,dimension(:,:)::Uo,Wo
    !velocity
    real(kind=dp),pointer,dimension(:,:)::Um,Wm,U,W
    
    
    !adv compute
    real(kind=dp),pointer,dimension(:,:)::radv,zadv,radvo,zadvo
    !diffusion
    real(kind=dp),pointer,dimension(:,:)::rdif,zdif,rdifo,zdifo
    !source
    real(kind=dp),pointer,dimension(:,:)::Rsource,Zsource,Uap,Rsourceo,Zsourceo
    real(kind=dp),pointer,dimension(:,:)::Rcsf,Zcsf
    real(kind=dp)::Rcsfavge,Zcsfavge

    !density 
    real(kind=dp)::Wdic,Wvis,Adic,Avis
    real(kind=dp),pointer,dimension(:,:)::F,FN,FO,Fnode
    real(kind=dp),pointer,dimension(:,:)::Denc,Visc
    !time
    real(kind=dp),save::dt
    real(kind=dp)::time,outtime,Cdmax,tend
    !pressure
    real(kind=dp),pointer,dimension(:,:)::Po,Pm,P
    
    !vof
    real(kind=dp),pointer,dimension(:,:)::Cfl,Nor,Noz,Coer,Coez,Coero,Coezo
    
    !solver
    real(kind=dp),pointer,dimension(:,:)::PAW,PAE,PAB,PAT,PAP
    real(kind=dp),pointer,dimension(:,:)::PRHS
    real(kind=dp),pointer,dimension(:,:)::NAME
    integer::number

    !csf
    real(kind=dp),pointer,dimension(:,:)::FAW,FAE,FAB,FAT,FAP
    real(kind=dp),pointer,dimension(:,:)::Curvc,Nnormalcr,Nnormalcz
    !check
    integer::Cyc,itermax
    real(kind=dp)::mindt

    real(kind=dp),pointer,dimension(:,:,:)::Hr,HZ

    real(kind=dp)::Agx,Agz

    !damping
    real(kind=dp)::amplitude,frequency


    end module par