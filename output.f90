subroutine output(index_file)
  use par
  implicit none
  real(kind=dp)::X,Z,Uout,Wout,HH
  real(kind=dp)::Pini,Pooin,force,Pforce,Visoforce
  real(kind=dp)::height,rm
  real(kind=dp)::masssum,massintx,massintz,centerx,centerz,centervel(2)
  character(len=100)::filename
  integer::IO,index_file,MWW,Lfile

  PARAMETER(MWW=6)
  character*1 var_name(mww)*1,index_name*4
  character*8 title
  data var_name/"X","Z","U","W","P","vof"/
  write(index_name,'(I4.4)') index_file
  title="zvof"//index_name

  IO=100+cyc
  outtime=time
  masssum=0.
  massintx=0.
  massintz=0.
  centerx=0.
  centerz=0.
  centervel(1)=0.
  centervel(2)=0.
  height=0.
  rm=0.
  write(*,*) index_file,'output tecplot files','*************************'
!center.data to save bubble's center and center's velocity
      open(10,position='Append',file='center.data')
      do i=RS,RE 
        do k=ZS,ZE 
          masssum=masssum+dr*dz*RC(i)*F(i,k)
          centerz=centerz+ZC(k)*dr*dz*RC(i)*F(i,k)
        end do 
      end do 
      centerz=centerz/masssum

      !height at the center/2
      i=RE/2
      height=0.
      do k=ZS,ZE 
        height=height+F(i,k)*dz
        if(F(i,k)>0..and.F(i,k)<1.0) then
          Pini=P(i,k)
        end if 
      end do 
      !force on the wall 
      i=RE 
      Pooin=0.
      force=0.
      Pforce=0.
      visoforce=0.
      do k=ZS,ZE
        Pforce=Pforce+(P(i,k)+Pooin)*RP(i)*dz*2.0*pi
        visoforce=visoforce+(Wvis*F(i,k)+Avis*(1.-F(i,k)))*(U(i+1,k)-U(i,k))/dr*RP(i)*dz*2.0*pi
        force=force+Pforce+Visoforce
      end do 

      write(10,100)time,centerz,height,Pini,force,Pforce,visoforce
100   format(7(f20.8))
      close(10)
      
      write(filename,"('file-'F8.4'.dat')")time

      open(unit=3,file=title//'.dat')
      write(3,*)"variables=",(',"',var_name(Lfile),'"',Lfile=1,MWW)
      write(3,110)NR-1,NZ-1,time
110   format('Zone i=',I5,',k=',I5,',F=POINT,solutiontime=',F8.4)
      do k=1,NZ-1
          do i=1,NR-1
              Uout=(U(i,k)+U(i+1,k))*0.5_dp
              Wout=(W(i,k)+W(i,k+1))*0.5_dp
              X=RC(i)
              Z=ZC(k)
              write(3,120)X,Z,Uout,Wout,P(i,k),F(i,k)
120           format(6(F20.8))
          end do
      end do
      close(3)
      !write the heigh data
      write(filename,"('surface elevatio file-'F8.4'.dat')")time
      open(unit=4,file=filename)
      write(4,*)'variables = "X","H"'
      write(4,130)NR-1,time
130   format('zone I=',I5,',F=POINT,solutiontime=',F8.4) 

      do i=1,NR-1
         X=RC(i)
         height=0.
         do k=1,NZ-1
            height=height+F(i,k)*dz
         end do 
         write(4,140)X,height
140      format(2(F20.8))  
      end do 
      close(4)

  return
  end subroutine output

subroutine save_data
    use par 
    implicit none
    open(20,position='rewind',file='save.data')
    write(20,*)nout
    write(20,*)time
    do k=1,NZ-1
            do i=1,NR
                write(20,*)U(i,k)
            end do    
        end do

    do i=1,NR-1
            do k=1,NZ 
                write(20,*)W(i,k)
            end do 
    end do 

    do i=1,NR-1
            do k=1,NZ-1
                write(20,*)P(i,k),F(i,k)
            end do    
    end do


    close(20)
    return 
    end subroutine save_data
