subroutine truncation(NWF,Na1,Na2,Na3,threshold,KS_R)
  implicit none 
  integer,intent(in)::NWF,Na1,Na2,Na3 
  real(8),intent(in)::threshold 
  complex(8),intent(inout)::KS_R(NWF,NWF,-Na1:Na1,-Na2:Na2,-Na3:Na3)   
  complex(8)::tij
  integer::ia1,ia2,ia3,ib,jb
  real(8),parameter::au=27.21151d0 
 !
 !OPEN(122,W,FILE='./dir-tr/dat.tr') 
 !
 OPEN(122,FILE='./dir-tr/dat.tr') 
 rewind(122)
  do ia1=-Na1,Na1
   do ia2=-Na2,Na2
    do ia3=-Na3,Na3
     do ib=1,NWF
      do jb=1,NWF 
       tij=KS_R(ib,jb,ia1,ia2,ia3)
       if(threshold>abs(tij))then 
        tij=0.0d0 
       else 
        write(122,'(5i5,2f15.7)') ib,jb,ia1,ia2,ia3,dble(tij)*au,abs(tij)*au
       endif 
       KS_R(ib,jb,ia1,ia2,ia3)=tij 
      enddo!jb 
     enddo!ib 
    enddo 
   enddo 
  enddo 
return
end
