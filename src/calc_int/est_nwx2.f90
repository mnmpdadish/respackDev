integer function algn235(inr)
implicit none
integer,intent(in):: inr
integer:: nr,m2,m3,m5,info
nr=inr
call fctck(nr,m2,m3,m5,info)
do while (info .eq. 1)
   nr = nr + 1
   call fctck(nr,m2,m3,m5,info)
end do
algn235 = nr 
return
end function algn235
!---
subroutine fctck(n,m2,m3,m5,info)
implicit none
integer,intent(in):: n
integer,intent(out):: m2, m3, m5, info
integer:: i
i=n
m2 = 0
m3 = 0
m5 = 0
info = 0
do while (i .ne. 1)
   if (mod(i,2) .eq. 0) then
      m2 = m2 + 1
      i = i / 2
   else if (mod(i,3) .eq. 0) then
      m3 = m3 + 1
      i = i / 3
   else if (mod(i,5) .eq. 0) then
      m5 = m5 + 1
      i = i / 5
   else
      info = 1
      exit
   end if
end do
return
end subroutine fctck
!---
