subroutine get_height_function_Hz
	use par 
    implicit none
    integer::ii
    real(kind=dp)::top1,down1,middle,top2,down2,top3,down3

    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    do i=RS,RE
            do k=ZS,ZE
            	do ii=-1,1
            		if(0.<F(i,k).and.F(i,k)<1.) then!in the transaint area
            			middle=F(i+ii,k)
            			top1=F(i+ii,k+1)
                        top2=F(i+ii,k+2)
                        top3=F(i+ii,k+3)
            			down1=F(i+ii,k-1)
                        down2=F(i+ii,k-2)
                        down3=F(i+ii,k-3)
            			Hz(i,k,ii)=(middle+top1+down1+top2+down2+top3+down3)*dz
            		end if 
            	end do !end do ii
            end do 
    end do 

    return
 end subroutine get_height_function_Hz


subroutine get_height_function_Hr
	use par 
    implicit none
    integer::jj,iter 
    real(kind=dp)::right1,middle,left1,right2,left2,right3,left3

    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    do i=RS,RE
            do k=ZS,ZE
            	do jj=1,-1,-1
            		if(0.<F(i,k).and.F(i,k)<1.) then
            			middle=F(i,k+jj)
            			right1=F(i+1,k+jj)
                        right2=F(i+2,k+jj)
                        right3=F(i+3,k+jj)
            			left1=F(i-1,k+jj)
                        left2=F(i-2,k+jj)
                        left3=F(i-3,k+jj)
            			Hr(i,k,jj)=(middle+right1+right2+right3+left1+left2+left3)*dr
            		end if 
            	end do 
            end do 
      end do 

      return
  end subroutine get_height_function_Hr


subroutine compute_curvature
    use par 
    implicit none 

    real(kind=dp)::nrx,nrz,HZrr,Hzr,Hrzz,Hrz 

    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    !compute curve
    do i=RS,RE 
        do k=ZS,ZE 
            if(0.<F(i,k).and.F(i,k)<1.) then
                nrx=(F(i+1,k)-F(i-1,k))/dr*0.5_dp
                nrz=(F(i,k+1)-F(i,k-1))/dz*0.5_dp
                if(abs(nrx)<abs(nrz)) then  !bootom to top H(r)
                    HZrr=(Hz(i,k,1)-2.0_dp*Hz(i,k,0)+Hz(i,k,-1))/dr**2
                    Hzr=(Hz(i,k,1)-Hz(i,k,-1))/dr*0.5_dp
                    Curv(i,k)=Hzrr/((sqrt(1+Hzr**2))**3+1.0e-20)+Hzr/(sqrt(1+Hzr**2)+1.0e-20)/RC(i)
                    Curv(i,k)=-Curv(i,k)
                    !if(Hzrr*nrz>=0) then
                        !Curv(i,k)=abs(Curv(i,k))
                    !else
                        !Curv(i,k)=-abs(Curv(i,k))
                    !end if 
                else !nrx>nrz right to left H(z)
                    Hrzz=(Hr(i,k,1)-2.0_dp*Hr(i,k,0)+Hr(i,k,-1))/dz**2
                    Hrz=(Hr(i,k,1)-Hr(i,k,-1))/dz*0.5_dp
                    Curv(i,k)=Hrzz/((sqrt(Hrz**2+1.0))**3+1.e-20)-1.0/Hr(i,k,0)/(sqrt(Hrz**2+1.0)+1.0e-20)
                    Curv(i,k)=-Curv(i,k)
                   ! if(Hrzz*nrx>=0) then
                        !Curv(i,k)=abs(Curv(i,k))
                    !else
                        !Curv(i,k)=abs(Curv(i,k))
                    !end if 
                end if 
            else
                Curv(i,k)=0.
            end if 
        end do
    end do

    return 

end subroutine compute_curvature

subroutine compute_csf
    use par 
    implicit none 

    real(kind=dp)::w1,W2,Cur,delta 

    RS=2
    RE=NR-1
    ZS=1
    ZE=NZ-1
    
    do i=RS,RE
            do k=ZS,ZE
                W1=F(i,k)*(1._dp-F(i,k))
                W2=F(i-1,k)*(1._dp-F(i-1,k))

                delta=(F(i,k)-F(i-1,k))/dr

                Cur=(W1*Curv(i,k)+W2*Curv(i-1,k))/(W1+W2+1.e-20_dp)

                Rcsf(i,k)=epiron*delta*Cur
            end do
    end do
    
    
    RS=1
    RE=NR-1
    
    ZS=2
    ZE=NZ-1
    
    do i=RS,RE
            do k=ZS,ZE
                W1=F(i,k)*(1._dp-F(i,k))
                W2=F(i,k-1)*(1._dp-F(i,k-1))

                delta=(F(i,k)-F(i,k-1))/dz

                Cur=(W1*Curv(i,k)+W2*Curv(i,k-1))/(W1+W2+1.e-20_dp)

                Zcsf(i,k)=epiron*delta*Cur
            end do
    end do

    return 
end subroutine compute_csf 

subroutine csf 
    use par
    implicit none
    call get_height_function_Hz
    call get_height_function_Hr
    call compute_curvature
    call compute_csf
    return
end subroutine csf 

