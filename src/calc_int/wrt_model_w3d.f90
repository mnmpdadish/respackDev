subroutine wrt_model(Na1,Na2,Na3,n_occ,UR,WR) 
implicit none 
integer::Na1,Na2,Na3,n_occ 
integer::N_element
integer::ia1,ia2,ia3,ib,jb,i 
integer,allocatable::unit_vec(:)!unit_vec(N_element) 
complex(8)::UR(n_occ,n_occ,-Na1:Na1,-Na2:Na2,-Na3:Na3)
real(8)::WR(-Na1:Na1,-Na2:Na2,-Na3:Na3)
real(8),parameter::au=27.21151d0
!
N_element=(2*Na1+1)*(2*Na2+1)*(2*Na3+1)  
!
!OPEN(301,W,FILE='zvo_ur.dat') 
!
OPEN(301,FILE='zvo_ur.dat') 
write(301,'(a)')'wannier90 format for vmcdry.out or HPHI -sdry'
write(301,'(i10)') n_occ
write(301,'(i10)') N_element 
!
allocate(unit_vec(N_element)); unit_vec=1
write(301,'(15i5)')(unit_vec(i),i=1,N_element) 
deallocate(unit_vec) 
do ia1=-Na1,Na1
 do ia2=-Na2,Na2
  do ia3=-Na3,Na3 
   do ib=1,n_occ
    do jb=1,n_occ
     write(301,'(i5,i5,i5,i5,i5,2f15.10)') ia1,ia2,ia3,ib,jb,UR(ib,jb,ia1,ia2,ia3)*WR(ia1,ia2,ia3)*au   
    enddo!jb
   enddo!ib
  enddo!ia3
 enddo!ia2
enddo!ia1 
return
end 
