  SUBROUTINE moving_average(width,length,band,mask,mask2,wsize,data)
    INTEGER,INTENT(IN)::width
    INTEGER,INTENT(IN)::length
    INTEGER,INTENT(IN)::band
    INTEGER,INTENT(IN)::wsize
    REAL(4)::dum(width,length,band)
    REAL(4),INTENT(INOUT) :: data(width,length,band)
    REAL(4),INTENT(IN) :: mask,mask2
    INTEGER::i,j,k,ii,jj,i3,j3,inum
    dum(:,:,:)=0.
    do j =1,length
       do i =1,width
!          if((data(i,j,1).ge.-0.1).and.(data(i,j,1).le.0.1)) then
!          print *,data(i,j,1)
!          endif
          if(data(i,j,1).ne.int(mask)) then
             inum=0
          do jj=1,wsize
             do ii=1,wsize
                i3=i+ii-(wsize-1)/2-1
                j3=j+jj-(wsize-1)/2-1
                if((i3.lt.1).or.(i3.gt.width)) cycle
                if((j3.lt.1).or.(j3.gt.length)) cycle
                if(data(i3,j3,1).eq.mask) cycle
                inum=1+inum            
                do k=1,band            
                   dum(i,j,k)=data(i3,j3,k)+dum(i,j,k)
                enddo
             enddo
          enddo
          dum(i,j,:)=dum(i,j,:)/real(inum)
!          dum(i,j,2)=dum(i,j,2)/real(inum)
!          dum(i,j,3)=dum(i,j,3)/real(inum)
       else
          dum(i,j,:)=mask2
          endif
       enddo
    enddo
    data(:,:,:)=dum(:,:,:)
  END SUBROUTINE moving_average
