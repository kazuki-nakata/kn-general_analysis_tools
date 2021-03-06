SUBROUTINE RESHAPE_BACKFORWARD_3D(data_i,data_o,ksize,jsize,isize)
    IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER(4),INTENT(IN) :: isize,jsize,ksize
  REAL(4),INTENT(IN) :: data_i(ksize,jsize,isize)
  REAL(4),INTENT(OUT) :: data_o(isize,jsize,ksize)
  do k=1,ksize
     do j=1,jsize
        do i =1,isize
     data_o(i,j,k)=data_i(k,j,i)
  enddo
enddo
enddo
END SUBROUTINE RESHAPE_BACKFORWARD_3D

  SUBROUTINE RESHAPE_BACKFORWARD_2D(data_i,data_o,jsize,isize)
    IMPLICIT NONE
  INTEGER :: i,j
  INTEGER(4),INTENT(IN) :: isize,jsize
  INTEGER(4),INTENT(IN) :: data_i(jsize,isize)
  INTEGER(4),INTENT(OUT) :: data_o(isize,jsize)

     do j=1,jsize
        do i =1,isize
           data_o(i,j)=real(data_i(j,i))
  enddo
enddo

END SUBROUTINE RESHAPE_BACKFORWARD_2D

  SUBROUTINE RESHAPE_BACKFORWARD_2D_REAL(data_i,data_o,jsize,isize)
    IMPLICIT NONE
  INTEGER :: i,j
  INTEGER(4),INTENT(IN) :: isize,jsize
  REAL(4),INTENT(IN) :: data_i(jsize,isize)
  REAL(4),INTENT(OUT) :: data_o(isize,jsize)

     do j=1,jsize
        do i =1,isize
           data_o(i,j)=real(data_i(j,i))
  enddo
enddo

END SUBROUTINE RESHAPE_BACKFORWARD_2D_REAL

SUBROUTINE RESHAPE_BACKFORWARD_4D(data_i,data_o,lsize,ksize,jsize,isize)
    IMPLICIT NONE
  INTEGER :: i,j,k,l
  INTEGER(4),INTENT(IN) :: isize,jsize,ksize,lsize
  REAL(4),INTENT(IN) :: data_i(lsize,ksize,jsize,isize)
  REAL(4),INTENT(OUT) :: data_o(isize,jsize,ksize,lsize)
  do l=1,lsize
  do k=1,ksize
     do j=1,jsize
        do i =1,isize
     data_o(i,j,k,l)=data_i(l,k,j,i)
  enddo
enddo
enddo
enddo

END SUBROUTINE RESHAPE_BACKFORWARD_4D
