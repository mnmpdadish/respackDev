subroutine wrt_model(Na1,Na2,Na3,n_occ,JR,WR) 
implicit none 
integer::Na1,Na2,Na3,n_occ 
integer::N_element 
integer::ia1,ia2,ia3,ib,jb,i 
integer,allocatable::unit_vec(:)!unit_vec(N_element) 
complex(8)::JR(n_occ,n_occ,-Na1:Na1,-Na2:Na2,-Na3:Na3)
real(8)::WR(-Na1:Na1,-Na2:Na2,-Na3:Na3)
real(8),parameter::au=27.21151d0
!
N_element=(2*Na1+1)*(2*Na2+1)*(2*Na3+1)  
!
!OPEN(302,W,FILE='zvo_jr.dat') 
!
OPEN(302,FILE='zvo_jr.dat') 
write(302,'(a)')'wannier90 format for vmcdry.out or HPHI -sdry'
write(302,'(i10)') n_occ
write(302,'(i10)') N_element 
!
allocate(unit_vec(N_element)); unit_vec=1
write(302,'(15i5)')(unit_vec(i),i=1,N_element) 
deallocate(unit_vec) 
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
