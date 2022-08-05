MODULE Calc
  IMPLICIT NONE
  INTEGER,PARAMETER :: null = 999
  REAL(4),PARAMETER :: Undef = 999
CONTAINS
  
SUBROUTINE Inner_Product(inval,isize2,jsize2,ksize2,inlat,inlon,isize,jsize,ymin,xmin,yint,xint,output)
   IMPLICIT NONE
   INTEGER :: i,j,ii,jj,k
   iNTEGER :: ipos,jpos,ipos2,jpos2
   REAL(4) :: dx,dy
   INTEGER(4),INTENT(IN) :: jsize,isize,isize2,jsize2,ksize2
   REAL(4),INTENT(IN) :: xmin,ymin,xint,yint
   REAL(4),INTENT(IN) :: inlat(1:isize,1:jsize),inlon(1:isize,1:jsize)
   REAL(4),INTENT(IN) :: inval(1:ksize2,1:isize2,1:jsize2)
   REAL(4) :: rlon(1:jsize2),rlat(1:isize2)
   REAL(4),INTENT(OUT):: output(1:ksize2,1:isize,1:jsize)

   rlon(1:jsize2) = (/ (xmin+xint*(i-1),i=1,jsize2) /)
   where(rlon > 360.)
      rlon=rlon-360.
   endwhere
   rlat(1:isize2) = (/ (ymin+yint*(i-1),i=1,isize2) /)
   print *,isize,jsize,isize2,jsize2,xmin,ymin,xint,yint
   
   do j = 1, jsize
      do i = 1, isize
         jpos=-999
         do jj = 2, jsize2
            if((inlon(i,j) > rlon(jj-1)) .and. &
                  (inlon(i,j) <= rlon(jj))) then
               jpos=jj
               jpos2=jj-1
            endif
         enddo


      enddo
END SUBROUTINE Inner_Product 

function bilin(dx,dy,Undef,var1,var2,var3,var4)
  implicit none
  REAL(4) :: bilin
  REAL(4),INTENT(in) :: dx, dy, Undef
  REAL(4),INTENT(in) :: var1,var2,var3,var4
!   REAL(4),INTENT(in) :: la1,la2,la3,la4

  if((var1 /= Undef) .and. (var2 /= Undef) .and. &
       & (var3 /= Undef) .and. (var4 /= Undef)) then
    bilin = (1.-dx)*(1.-dy)*var1 &
         & +    dx *(1.-dy)*var2 &
         & +    dx *    dy *var3 &
         & +(1.-dx)*    dy *var4
  else
    bilin = Undef
  endif
end function bilin

END MODULE Calc
