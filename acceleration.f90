subroutine acceleration(nowtime,aga,agb)
    use par
    implicit none
    real(kind=dp),intent(in)::nowtime
    real(kind=dp),intent(out)::aga,agb
    aga=0._dp
    agb=-9.81*amplitude*abs(sin(2.0*pi*frequency*time))
!     agb=-9.81*(1.0+amplitude*sin(2.0*pi*frequency*time))
  
    !agb=0._dp
    return
    end subroutine acceleration
    
   real function  sgn(numbera)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp)::numbera
    if(numbera>0._dp) then
        sgn=1.
    else if(numbera==0._dp) then
        sgn=0.
    else
        sgn=-1.
    end if
    end function sgn