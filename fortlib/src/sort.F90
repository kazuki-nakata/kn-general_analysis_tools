subroutine heapsort(n,array,pairs1,pairs2)
  implicit none
  integer(4),intent(in) :: n
  real(4),intent(inout) :: array(1:n)
  REAL(4),intent(inout) :: pairs1(1:n),pairs2(1:n)
  REAL(4) :: dum_pairs1,dum_pairs2
  integer ::i,k,j,l
!  real(4) :: t
  double precision :: t
 
  if(n.le.0)then
     write(6,*)"Error, at heapsort"; stop
  endif
  if(n.eq.1)return

  l=n/2+1
  k=n
  do while(k.ne.1)
     if(l.gt.1)then
        l=l-1
        t=array(L)
        dum_pairs1=pairs1(L)
        dum_pairs2=pairs2(L)
     else
        t=array(k)
        dum_pairs1=pairs1(k)
        dum_pairs2=pairs2(k)
        
        array(k)=array(1)
        pairs1(k)=pairs1(1)
        pairs2(k)=pairs2(1)        
        k=k-1
        if(k.eq.1) then
           array(1)=t
           pairs1(1)=dum_pairs1
           pairs2(1)=dum_pairs2
           exit
        endif
     endif
     i=l
     j=l+l
     do while(j.le.k)
        if(j.lt.k)then
           if(array(j).lt.array(j+1))j=j+1
        endif
        if (t.lt.array(j))then
           array(i)=array(j)
           pairs1(i)=pairs1(j)
           pairs2(i)=pairs2(j)             
           i=j
           j=j+j
        else
           j=k+1
        endif
     enddo
     array(i)=t
     pairs1(i)=dum_pairs1
     pairs2(i)=dum_pairs2     
  enddo

  return
end subroutine heapsort
