subroutine vof_explicit
	use par
	implicit none
	call vof_cfl
	call vof_no
	call vof_co
	call vof_update
	return
end subroutine vof_explicit

subroutine vof_cfl   !cf
    use par
    implicit none
    
    real(kind=dp)::vp
    RS=1
    RE=NR-1
    ZS=1
    ZE=NZ-1
    !....................................................................
    do i=RS,RE
            do k=ZS,ZE
                if(cylinder==1) then
                    vp=RC(i)*dr*dz
                     Cfl(i,k)=max(-U(i,k)*RP(i)*dz/vp*dt,0._dp)+max(U(i+1,k)*RP(i+1)*dz/vp*dt,0._dp)+&
                        max(-W(i,k)*RC(i)*dr/vp*dt,0._dp)+max(W(i,k+1)*RC(i)*dr/vp*dt,0._dp)
                else 
                    vp=dr*dz
                    Cfl(i,k)=max(-U(i,k)*dz/vp*dt,0._dp)+max(U(i+1,k)*dz/vp*dt,0._dp)+&
                        max(-W(i,k)*dr/vp*dt,0._dp)+max(W(i,k+1)*dr/vp*dt,0._dp)
                end if 
            end do
    end do
    return
end subroutine vof_cfl

subroutine vof_no   !F½ü±ÚÃæ0ÌÝ¶È¼ÙÉè
    use par
    implicit none
    real(kind=dp)::Fwf,Fef,Fbf,Ftf
    real(kind=dp)::Awf,Aef,Abf,Atf
    real(kind=dp),pointer::DfDr(:,:),DfDz(:,:)
    real(kind=dp)::vp,NN
    
    allocate(DfDr(1:NR,1:NZ),DfDz(1:NR,1:NZ))
    do i=1,NR
            do k=1,NZ
                DfDr(i,k)=0.
                DfDz(i,k)=0.
            end do
    end do
    
    
    
    RS=1
    RE=NR-1
    ZS=1
    ZE=NZ-1

    !Fn first update by csf
    !then update by vof_update
    do i=1,NR 
        do k=1,NZ 
            Fnode(i,k)=0.25*(FN(i,k)+FN(i-1,k)+FN(i,k-1)+FN(i-1,k-1))
        end do 
    end do 
    
    do i=RS,RE
            do k=ZS,ZE
                if(cylinder==1) then 
                    vp=RC(i)*dr*dz
                    Awf=RP(i)*dz
                    Aef=RP(i+1)*dz
                    Abf=RC(i)*dr
                    Atf=RC(i)*dr
                else 
                    vp=dr*dz
                    Awf=dz
                    Aef=dz
                
                    Abf=dr
                    Atf=dr
                end if 

                Fwf=0.5_dp*(Fnode(i,k)+Fnode(i,k+1))
                Fef=0.5_dp*(Fnode(i+1,k)+Fnode(i+1,k+1))
                Fbf=0.5_dp*(Fnode(i,k)+Fnode(i+1,k))
                Ftf=0.5_dp*(Fnode(i,k+1)+Fnode(i+1,k+1))
                DfDr(i,k)=(Fef-Fwf)/dr
                DfDz(i,k)=(Ftf-Fbf)/dz

                NN=sqrt(DfDr(i,k)*DfDr(i,k)+Dfdz(i,k)*Dfdz(i,k))+1.e-10_dp
                Nor(i,k)=DfDr(i,k)/NN
                Noz(i,k)=Dfdz(i,k)/NN
            end do
    end do
                
    !boundary   Nor,Nos£¬Noz in these additional cells will be use
     do k=ZS,ZE
            i=RS-1
            if(cylinder==1) then
                Nor(i,k)=-Nor(RS,k)
            else
                Nor(i,k)=Nor(RS,k)
            end if 
            i=RE+1
            Nor(i,k)=Nor(RE,k)
    end do
    
    do i=RS-1,RE+1
            k=ZS-1
            Noz(i,k)=Noz(i,ZS)!Noz(i,ZS)
            k=ZE+1
            Noz(i,k)=Noz(i,ZE)
    end do
    
    
    deallocate(DfDr,DfDz)
    
    return
    end  subroutine vof_no

subroutine vof_co
    use par
    implicit none
    real(kind=dp)::Uf,Fi(4),No(2),Cf(2),Fcoe
    RS=1
    RE=NR-1

    ZS=1
    ZE=NZ-1
    do k=ZS,ZE  !i=1ÒÔ¼°i=NR´¦µÄÃæ¿ÉÒÔ²»Ëã
            do i=2,RE
                Uf=U(i,k)
                Fi(1)=F(i-2,k)
                Fi(2)=F(i-1,k)
                Fi(3)=F(i,k)
                Fi(4)=F(i+1,k)
                No(1)=Nor(i-1,k)
                No(2)=Nor(i,k)
                Cf(1)=Cfl(i-1,k)
                Cf(2)=Cfl(i,k)
                call vof_HR(Uf,Fi,Cf,No,Fcoe)
                Coer(i,k)=Fcoe
            end do
    end do
    
    
    do i=RS,RE
            do k=2,ZE !k==1,and k==NZ
                Uf=U(i,k)
                Fi(1)=F(i,k-2)
                Fi(2)=F(i,k-1)
                Fi(3)=F(i,k)
                Fi(4)=F(i,k+1)
                No(1)=Noz(i,k-1)
                No(2)=Noz(i,k)
                Cf(1)=Cfl(i,k-1)
                Cf(2)=Cfl(i,k)
                call vof_HR(Uf,Fi,Cf,No,Fcoe)
                Coez(i,k)=Fcoe
            end do
    end do
    
    return
end subroutine vof_co

subroutine vof_update  
    use par 
    implicit none
    integer::vofstep,it,itin,iter,itmax
    integer::bound,boundsgn
    real(kind=dp)::vel_vol_correct
    real(kind=dp)::err,vp
    real(kind=dp)::Ff,Ffo
    real(kind=dp)::FD,FA,FDO,FAO
    real(kind=dp)::totalvol
    real(kind=dp)::divv,dtt
    real(kind=dp)::Phi(4),face_vlome(4),face_vlomeout(4),volumesuper
    real(kind=dp),pointer,dimension(:,:)::fluxr,fluxz,vof_face_volr,vof_face_volz

    real(kind=dp),pointer,dimension(:,:)::Tempnew 

    allocate(fluxr(1:NR,1:NZ),fluxz(1:NR,1:NZ))
    allocate(tempnew(1:NR,1:NZ))
    allocate(vof_face_volr(1:NR,1:NZ),vof_face_volz(1:NR,1:NZ))


    vofstep=5
    itmax=100
    bound=1
    dtt=0.67_dp*dt
    RS=1
    RE=NR-1
    
    ZS=1
    ZE=NZ-1
    do i=RS,RE
            do k=ZS,ZE
                Coero(i,k)=Coer(i,k)
                Coezo(i,k)=Coez(i,k)
            end do
    end do

    !.........................dual time step
    it=1
    err=1.0_dp
    do while(it<vofstep.and.err>1.e-6)
        !.........................
        it=it+1
        err=0._dp
        !update the face volume
        do i=1,NR 
            do k=1,NZ-1
                if(cylinder==1) then
                    fluxr(i,k)=U(i,k)*RP(i)*dz
                else
                    fluxr(i,k)=U(i,k)*dz
                end if 
                if(fluxr(i,k)>0.) then
                    FD=F(i-1,k)
                    FA=F(i,k)
                else
                    FA=F(i-1,k)
                    FD=F(i,k)
                end if 
                if(fluxr(i,k)>0.) then
                    FDO=Fo(i-1,k)
                    FAO=Fo(i,k)
                else
                    FAO=Fo(i-1,k)
                    FDO=Fo(i,k)
                end if 
                if(i==RS.or.i==NR) then
                    Ff=0.5_dp*(F(i,k)+F(i-1,k))
                    Ffo=0.5_dp*(Fo(i,k)+Fo(i-1,k))
                else
                    Ff=(1._dp-Coer(i,k))*FD+Coer(i,k)*FA
                    Ffo=(1._dp-Coero(i,k))*FDO+Coero(i,k)*FAO
                end if
                vof_face_volr(i,k)=dtt*0.5_dp*(Ff*fluxr(i,k)+Ffo*fluxr(i,k))
            end do 
        end do 

        do k=1,NZ
            do i=1,NR-1
                if(cylinder==1) then
                    fluxz(i,k)=W(i,k)*RC(i)*dr 
                else
                    fluxz(i,k)=W(i,k)*dr 
                end if 
                if(fluxz(i,k)>0.) then
                    FD=F(i,k-1)
                    FA=F(i,k)
                else
                    FD=F(i,k)
                    FA=F(i,k-1)
                end if
                if(fluxz(i,k)>0.) then
                    FDO=FO(i,k-1)
                    FAO=FO(i,k)
                else
                    FDO=FO(i,k)
                    FAO=FO(i,k-1)
                end if

                if(k==ZS.or.k==NZ) then
                    Ff=0.5_dp*(F(i,k)+F(i,k-1))
                    Ffo=0.5_dp*(Fo(i,k)+Fo(i,k-1))
                else
                    Ff=(1._dp-Coez(i,k))*FD+Coez(i,k)*FA
                    Ffo=(1._dp-Coezo(i,k))*FDO+Coezo(i,k)*FAO
                end if
                vof_face_volz(i,k)=dtt*0.5_dp*(Ff*fluxz(i,k)+Ffo*fluxz(i,k))
            end do 
        end do 

        !setp 2 start to update vof values
        RS=1
        RE=NR-1
    
        ZS=1
        ZE=NZ-1
        do i=RS,RE 
            do k=ZS,ZE 
                if(cylinder==1) then
                    Vp=dr*dz*Rc(i)
                else
                    Vp=dr*dz
                end if 
                totalvol=vof_face_volr(i+1,k)-vof_face_volr(i,k)+vof_face_volz(i,k+1)-vof_face_volz(i,k)
                Tempnew(i,k)=F(i,k)-totalvol/vp-(F(i,k)-FO(i,k))/dt*dtt

                vel_vol_correct=(fluxr(i+1,k)-fluxr(i,k)+fluxz(i,k+1)-fluxz(i,k))*dtt/vp
                if(tempnew(i,k)>1.0) then
                    tempnew(i,k)=tempnew(i,k)-abs(vel_vol_correct)
                else if(tempnew(i,k)<1.0) then
                    tempnew(i,k)=tempnew(i,k)+abs(vel_vol_correct)
                end if 
!                 TEMP=MIN(1.,MAX(TEMP,0.))
                if(tempnew(i,k)<0..or.tempnew(i,k)>1.) then
                    boundsgn=0
                else
                    boundsgn=1
                end if 
                bound=bound*boundsgn
            end do 
        end do 
        !check bounded and correxted it 
        !now get tempnew(i,k),tempold=F(i,k),FO(i,k)
        call boundary_f 
        iter=1
        do while(bound==0.and.iter<=itmax)
            !tempF is false 
!             write(*,*)'corrcting',iter 
            bound=1
            iter=iter+1
            RS=1
            RE=NR-1
    
            ZS=1
            ZE=NZ-1
            do i=RS,RE
                do k=ZS,ZE
                    if(tempnew(i,k)<0..or.tempnew(i,k)>1.) then
                        if(cylinder==1) then
                            vp=Rc(i)*dr*dz
                        else
                            vp=dr*dz 
                        end if 
                        do itin=1,2
                            Phi(itin)=fluxr(i+itin-1,k)
                            Phi(itin+2)=fluxz(i,k+itin-1)
                            face_vlome(itin)=vof_face_volr(i+itin-1,k)
                            face_vlome(itin+2)=vof_face_volz(i,k+itin-1)
                        end do 
                        if(tempnew(i,k)<0.) then
                            volumesuper=tempnew(i,k)*Vp
                        else if(tempnew(i,k)>1.) then
                            volumesuper=(tempnew(i,k)-1.)*Vp
                        end if 
                        call face_volume_correct(Phi,face_vlome,volumesuper,dtt,face_vlomeout)
                        do itin=1,2
                            vof_face_volr(i+itin-1,k)=face_vlomeout(itin)
                            vof_face_volz(i,k+itin-1)=face_vlomeout(itin+2)
                        end do 
                    end if 
                end do 
            end do 
            !recompute the vof_values
            RS=1
            RE=NR-1
    
            ZS=1
            ZE=NZ-1
            do i=RS,RE 
                do k=ZS,ZE 
                    if(cylinder==1) then
                        Vp=dr*dz*Rc(i)
                    else
                        Vp=dr*dz
                    end if 
                    totalvol=vof_face_volr(i+1,k)-vof_face_volr(i,k)+vof_face_volz(i,k+1)-vof_face_volz(i,k)
                    Tempnew(i,k)=F(i,k)-totalvol/vp-(F(i,k)-FO(i,k))/dt*dtt

                    vel_vol_correct=(fluxr(i+1,k)-fluxr(i,k)+fluxz(i,k+1)-fluxz(i,k))*dtt/vp
                    if(tempnew(i,k)>1.0) then
                        tempnew(i,k)=tempnew(i,k)-abs(vel_vol_correct)
                    else if(tempnew(i,k)<1.0) then
                        tempnew(i,k)=tempnew(i,k)+abs(vel_vol_correct)
                    end if 
!                   TEMP=MIN(1.,MAX(TEMP,0.))
                    if(tempnew(i,k)<0..or.tempnew(i,k)>1.) then
                        boundsgn=0
                    else
                        boundsgn=1
                    end if
                    bound=bound*boundsgn
                end do 
            end do 

        end do !end of corect


        call boundary_f 
        !end of a update in one daul time 
        !check the err of temp and F(i,k)
        do i=RS,RE 
            do k=ZS,ZE 
                err=max(err,abs(Tempnew(i,k)-F(i,k)))
            end do 
        end do 
        !no matter how is the result F(i,k) will be updated
        do i=RS,RE 
            do k=ZS,ZE 
                Tempnew(i,k)=F(i,k)
            end do 
        end do 
        !prepare for the next n^ compute and NA^compute
        do i=RS,RE 
            do k=ZS,ZE
                 if(F(i,k)<0..or.F(i,k)>1.) then
                     write(*,*)'correctfalse'
                     write(*,*)i,k,F(i,k)
                 end if 
            end do 
        end do 
      
        do i=-1,NR+1  
                do k=-1,NZ+1
                    FN(i,k)=F(i,k)
                end do
        end do
        call smooth_f
        !end of prepare
        call vof_no
        call vof_co

    end do !end of update 
    deallocate(fluxr,fluxz,vof_face_volr,vof_face_volz,Tempnew )
    return
    
end subroutine vof_update

subroutine vof_HR(Ufin,Fiin,Cfin,Noin,Fcoeout)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::Ufin,Fiin(4),Noin(2),Cfin(2)
    real(kind=dp),intent(out)::Fcoeout
    real(kind=dp)::FUin,FDin,FAin,anglein,courantin
    real(kind=dp)::FDNin,FCBCin,FUQin
    if(Ufin>0._dp) then
        FUin=Fiin(1)
        FDin=Fiin(2)
        FAin=Fiin(3)
        anglein=Noin(1)
        courantin=Cfin(1)
    else
        FUin=Fiin(4)
        FDin=Fiin(3)
        FAin=Fiin(2)
        anglein=Noin(2)
        courantin=Cfin(2)
    end if
    
    if(abs(FAin-FUin)>1.E-20) then
        FDNin=(FDin-FUin)/(FAin-FUin)
    else
        FDNin=0._dp !´ËÊ±FfÈ¡FU£¬FD£¬FAÃ»ÓÐÇø±ð
    end if
    !CBC value
    if(FDNin>0._dp.and.FDNin<1._dp.and.abs(courantin)>1.e-20_dp) then
        FCBCin=min(FDNin/(courantin+1.e-30_dp),1._dp)
    else
        FCBCin=FDNin
    end if
    !UQ value
    if(FDNin>0._dp.and.FDNin<1._dp) then
        FUQin=min((8._dp*courantin*FDNin+(1._dp-courantin)*(6._dp*FDNin+3._dp))/8._dp,FCBCin)
    else
        FUQin=FDNin
    end if
    anglein=abs(anglein)
    anglein=acos(anglein)
    anglein=min((cos(2._dp*anglein)+1._dp)*0.5_dp,1._dp)
    Fcoeout=anglein*FCBCin+(1._dp-anglein)*FUQin
    
    if(abs(1._dp-FDNin)>1.e-20_dp) then
        Fcoeout=(Fcoeout-FDNin)/(1._dp-FDNin)
    else
        Fcoeout=0._dp
    end if
    
    Fcoeout=min(max(Fcoeout,0._dp),1._dp)
    return
    end subroutine vof_HR
                         
    subroutine smooth_f
    use par
    implicit none
    integer::it,itm
    real(kind=dp)::Ffe,Ffw,Ffb,Fft
    real(kind=dp)::Afe,Afw,Afb,Aft
    real(kind=dp),pointer,dimension(:,:)::Ftem
    allocate(Ftem(1:NR-1,1:NZ-1))
    
    it=1
    itm=3
    
    do i=1,NR-1
            do k=1,NZ-1
                Ftem(i,k)=0_dp
            end do
    end do

    !......................................
    do while(it<itm)
        it=it+1 
        RS=1
        RE=NR-1
        ZS=1
        ZE=NZ-1 
        do i=RS,RE
                do k=ZS,ZE
                    !the i,k respent the cell number
                    !get the smooth F
                    Ffw=0.5_dp*(FN(i,k)+FN(i-1,k))
                    Ffe=0.5_dp*(FN(i,k)+FN(i+1,k))
                    Ffb=0.5_dp*(FN(i,k)+FN(i,k-1))
                    Fft=0.5_dp*(FN(i,k)+FN(i,k+1))
                    if(cylinder==1) then 
                        Afw=dz*RP(i)
                        Afe=dz*RP(i+1)
                        Afb=dr*RC(i)
                        Aft=dr*RC(i)
                    else 
                        Afw=dz
                        Afe=dz
                        Afb=dr
                        Aft=dr
                    end if 
                    Ftem(i,k)=(Ffe*Afe+Ffw*Afw+Ffb*Afb+Fft*Aft)
                    Ftem(i,k)=Ftem(i,k)/(Afe+Afw+Afb+Aft)
                end do
        end do
        
        RS=1
        RE=NR-1
        ZS=1
        ZE=NZ-1 
        do i=RS,RE
                do k=ZS,ZE      
                    FN(i,k)=Ftem(i,k)
                end do
        end do
        
    end do
    
    !............................................................................ 
    !boundary condition use to compute the additional csf 
    call boundary_FN
    
    deallocate(Ftem)
    
    return
    end subroutine smooth_f


subroutine face_volume_correct(Phiin,F1volume,volumeplus,deltat,F1volumeout)
    !correct the volume pass the face during [t,t+\delta t]
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp)::Phiin(4),F1volume(4)
    logical::isdownwind(4)
    real(kind=dp)::volumeplus,deltat,testvol
    real(kind=dp),intent(out)::F1volumeout(4)
    real(kind=dp)::total,deliveredvolume,a
    integer::deliversize,it,i 
    logical::capacity(4)
    !transport to positive values
    Phiin(1)=-Phiin(1)
    F1volume(1)=-F1volume(1)
    Phiin(3)=-Phiin(3)
    F1volume(3)=-F1volume(3)
    !decide the values of ifisdownwind
    do i=1,4
        if(Phiin(i)>0.) then
            isdownwind(i)=.true.
        else
            isdownwind(i)=.false.
        end if 
    end do 
    !start
    total=0.
    deliversize=0
    do i=1,4
        if(isdownwind(i)) then
            capacity(i)=.true.
            total=total+Phiin(i)
            deliversize=deliversize+1
        else
            capacity(i)=.false.
        end if 
    end do 


    if(volumeplus>0.) then
        !start deliver
        it=1
        do while(volumeplus>0..and.it<=deliversize) 
            !volumeplus==0. deliver over
            !it==deliversize,all downwindface get its capacity
            deliveredvolume=0.
            testvol=0.
            do i=1,4
                if(capacity(i)) then
                    testvol=testvol+Phiin(i)*deltat-F1volume(i)
                end if 
            end do 
!             if(testvol<volumeplus) then
!                 write(*,*)testvol,volumeplus
!             end if 
            !start a new deliver with (vloumeplus,capacity,total)
            do i=1,4
                if(capacity(i)) then
                    a=min(Phiin(i)*deltat-F1volume(i),volumeplus*Phiin(i)/total)
                    deliveredvolume=deliveredvolume+a
                    F1volume(i)=F1volume(i)+a
                    if(Phiin(i)*deltat-F1volume(i)<=volumeplus*Phiin(i)/total) then 
                        !stitl use the flux as weight
                        total=total-Phiin(i)
                        capacity(i)=.false. !capacity reached
                    end if 
                end if 
            end do 
            !the rest volumeplus
            volumeplus=volumeplus-deliveredvolume
            it=it+1
        !end deliver
!             write(*,*)'volumeplus,it '
!             write(*,*)volumeplus,it 
        end do 
    else if(volumeplus<0.) then
        it=1
        do while(volumeplus<0..and.it<=deliversize)
        !start with (vloumeplus)
        deliveredvolume=0.!delivered in one cycle
            do i=1,4
                if(capacity(i)) then
                    a=max(-F1volume(i),volumeplus*Phiin(i)/total)
                    deliveredvolume=deliveredvolume+a
                    F1volume(i)=F1volume(i)+a
                    if(-F1volume(i)>=volumeplus*Phiin(i)/total)   then!reach the bottom
                        total=total-Phiin(i)
                        capacity(i)=.false.
                    end if 
                end if 
            end do 
            volumeplus=volumeplus-deliveredvolume
            it=it+1
        end do 
    end if 

  !convert to the original sign
  F1volume(1)=-F1volume(1)
  F1volume(3)=-F1volume(3)
  !convert over 
  do i=1,4
        F1volumeout(i)=F1volume(i)
  end do 

  return
end subroutine  face_volume_correct


