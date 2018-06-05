subroutine wrt_mvmc(NTK,Na1,Na2,Na3,n_occ,JR,WR) 
implicit none 
integer::NTK,Na1,Na2,Na3,n_occ 
integer::ia1,ia2,ia3,ib,jb,i 
integer::unit_vec(NTK) 
complex(8)::JR(n_occ,n_occ,-Na1:Na1,-Na2:Na2,-Na3:Na3)
real(8)::WR(-Na1:Na1,-Na2:Na2,-Na3:Na3)
real(8),parameter::au=27.21151d0
!
!OPEN(302,W,FILE='zvo_jr.dat') 
!
OPEN(302,FILE='zvo_jr.dat') 
write(302,'(a)')'wannier90 format for mvmcdry'
write(302,'(i10)') n_occ
write(302,'(i10)') NTK 
unit_vec=1
write(302,'(15i5)')(unit_vec(i),i=1,NTK) 
do ia1=-Na1,Na1
 do ia2=-Na2,Na2 
  do ia3=-Na3,Na3 
   do ib=1,n_occ
    do jb=1,n_occ
     write(302,'(i5,i5,i5,i5,i5,2f15.10)') ia1,ia2,ia3,ib,jb,JR(ib,jb,ia1,ia2,ia3)*WR(ia1,ia2,ia3)*au   
    enddo!jb
   enddo!ib
  enddo!ia3
 enddo!ia2
enddo!ia1 
return
end 
