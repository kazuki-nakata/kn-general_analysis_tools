subroutine heapsort2(n,array)
  implicit none
  integer(4),intent(in) :: n
  real(4),intent(inout) :: array(1:n)
  integer(4) :: dum_pairs1,dum_pairs2
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
     else
        t=array(k)
        
        array(k)=array(1)    
        k=k-1
        if(k.eq.1) then
           array(1)=t
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
           i=j
           j=j+j
        else
           j=k+1
        endif
     enddo
     array(i)=t
  enddo

  return
end subroutine heapsort2
