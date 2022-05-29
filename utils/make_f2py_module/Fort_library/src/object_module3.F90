MODULE OBJECT_MODULE3
  IMPLICIT NONE
  SAVE
!  PRIVATE
  !  PUBLIC init_unionfind,find_uf,union_uf,get_ipos_uf
  INTEGER(4),PARAMETER :: n=1
    TYPE djsdata_struct      !disjoint-set data   
     INTEGER(4) :: par
     INTEGER(4) :: size
     INTEGER(4) :: ipos
     INTEGER(4) :: jpos
     INTEGER(4)    :: property(n)
  END TYPE djsdata_struct
  
  TYPE  doubly_linked_list
     INTEGER(4) :: par               !order from a chi to the parent
     INTEGER(4) :: chi               !order from a parent to the child
  END TYPE doubly_linked_list
  
  TYPE  singly_linked_list
     INTEGER(4) :: dchi               !distant child node
     INTEGER(4) :: nchi               !neibor child node
  END TYPE singly_linked_list

  INTEGER(4) :: smax_uf
  
  INTEGER(4),ALLOCATABLE :: ptrt(:,:)
  TYPE(djsdata_struct), ALLOCATABLE :: djsdata(:)
     TYPE(doubly_linked_list), ALLOCATABLE :: dllptr(:)
     TYPE(singly_linked_list), ALLOCATABLE :: sllptr(:)     
CONTAINS

  SUBROUTINE init_unionfind(width,length,data,n_prop,mask,num)
    INTEGER,INTENT(IN)::width
    INTEGER,INTENT(IN)::length
    INTEGER,INTENT(IN)::n_prop
    INTEGER,INTENT(IN) :: data(width,length,n_prop)
    INTEGER,INTENT(IN) :: mask
    INTEGER,INTENT(OUT):: num
    INTEGER :: i, j, k
    INTEGER(4) :: ptr
    INTEGER(4) :: count
    

!----------------data count---------------------- 
    ptr=0

    do j =1,length
       do i =1,width
          if(data(i,j,1).eq.int(mask)) cycle
          ptr=1+ptr
       enddo
    enddo
    smax_uf=ptr
!------------------------------------------------
    
    ALLOCATE(ptrt(1:width,1:length))          
    ALLOCATE(djsdata(1:smax_uf))
    ALLOCATE(sllptr(1:smax_uf))  
    ALLOCATE(dllptr(1:smax_uf))

   ptrt(:,:)=0          
   ptr=0          
    do j =1,length
       do i =1,width
   if(data(i,j,1).eq.int(mask)) cycle

      ptr=ptr+1
      ptrt(i,j)=ptr
      djsdata(ptr)%par = ptr
    do k=1,n_prop
      djsdata(ptr)%property(k)=data(i,j,k)
    enddo
    djsdata(ptr)%ipos = i
    djsdata(ptr)%jpos = j
    djsdata(ptr)%size = 1
    dllptr(ptr)%chi = ptr-1
    dllptr(ptr)%par = ptr+1
    sllptr(ptr)%nchi = ptr
    sllptr(ptr)%dchi = ptr   

 enddo
enddo
smax_uf=ptr
num=smax_uf
  END SUBROUTINE init_unionfind
  
  SUBROUTINE union_uf(tmp1,tmp2,index_union)
    INTEGER(4) :: tmp1,tmp2
    INTEGER(4) :: ph1,ph2
    INTEGER(4) :: s1,s2
    INTEGER(4) :: i1,i2,i3,i4
    INTEGER(4) :: chi,par,dchi1,dchi2
    INTEGER(4),INTENT(OUT)::index_union
    ph1=find_uf(tmp1)
    ph2=find_uf(tmp2)
    !    print *,ph1,ph2
    if(ph1.eq.ph2) then
       index_union=999
       return !
    endif
    s1=djsdata(ph1)%size
    s2=djsdata(ph2)%size
    if(s1 .lt. s2) then
       i1=ph2;i2=ph1;i3=tmp2;i4=tmp1
    else
       i1=ph1;i2=ph2;i3=tmp1;i4=tmp2
    endif

    call set_par_uf(i2,i1)
   djsdata(i1)%size=s1+s2
   djsdata(i1)%property(1)=(s1*djsdata(ph1)%property(1)+s2*djsdata(ph2)%property(1))/(s1+s2)
   djsdata(i2)%size=0
   
   dchi1=get_dchi_sll(i1)
   dchi2=get_dchi_sll(i2)   
   call set_nchi_sll(dchi1,i2)
   call set_dchi_sll(i1, dchi2)
   call set_nchi_sll(dchi2,0)
   
   chi=get_chi_dll(i2)
   par=get_par_dll(i2)  
   if(chi.eq.0) then
      call set_chi_dll(par,chi)
!      call set_par_dll(i1,-1)      
   elseif(par.eq.smax_uf+1) then
      call set_par_dll(chi,par)
!      call set_par_dll(i2,-1)
   else
      call set_chi_dll(par,chi)
      call set_par_dll(chi,par)
   endif
!   print *,par,chi,ph2

    index_union=i1
!   print *,get_chi_dll(chi),get_par_dll(chi),'chi'
!   print *,get_chi_dll(par),get_par_dll(par),'par'   
  END SUBROUTINE union_uf

  INTEGER(4) FUNCTION get_size_uf(ph)
    INTEGER(4), INTENT(IN) :: ph
    get_size_uf = djsdata(ph)%size
  END FUNCTION get_size_uf
  
  INTEGER(4) FUNCTION get_ipos_uf(ph)
    INTEGER(4), INTENT(IN) :: ph
    get_ipos_uf = djsdata(ph)%ipos
  END FUNCTION get_ipos_uf
  
  REAL(4) FUNCTION get_property_uf(ph,vn)
    INTEGER(4), INTENT(IN) :: ph
    INTEGER(4), INTENT(IN) :: vn
    get_property_uf = djsdata(ph)%property(vn)
  END FUNCTION get_property_uf
  
  INTEGER(4) FUNCTION get_jpos_uf(ph)
    INTEGER(4), INTENT(IN) :: ph
    get_jpos_uf = djsdata(ph)%jpos
  END FUNCTION get_jpos_uf
  
  SUBROUTINE set_par_uf(ph1, ph2)
    INTEGER(4), INTENT(IN) :: ph1,ph2
    djsdata(ph1)%par = ph2
  END SUBROUTINE set_par_uf
  
  SUBROUTINE set_property_uf(ph1, vn, val)
    INTEGER(4), INTENT(IN) :: ph1
    INTEGER(4), INTENT(IN) :: vn
    INTEGER(4), INTENT(IN) :: val
    djsdata(ph1)%property(vn) = val
  END SUBROUTINE set_property_uf
!-------------------------------------------------------
  INTEGER(4) FUNCTION get_nchi_sll(ph)
    INTEGER(4), INTENT(IN) :: ph
    get_nchi_sll = sllptr(ph)%nchi
  END FUNCTION get_nchi_sll

  INTEGER(4) FUNCTION get_dchi_sll(ph)
    INTEGER(4), INTENT(IN) :: ph
    get_dchi_sll = sllptr(ph)%dchi
  END FUNCTION get_dchi_sll

  INTEGER(4) FUNCTION get_chi_dll(ph)
    INTEGER(4), INTENT(IN) :: ph
    get_chi_dll = dllptr(ph)%chi
  END FUNCTION get_chi_dll

  INTEGER(4) FUNCTION get_par_dll(ph)
    INTEGER(4), INTENT(IN) :: ph
    get_par_dll = dllptr(ph)%par
  END FUNCTION get_par_dll
    
  SUBROUTINE set_par_dll(ph1, ph2)
    INTEGER(4), INTENT(IN) :: ph1,ph2
    dllptr(ph1)%par = ph2
  END SUBROUTINE set_par_dll

  SUBROUTINE set_chi_dll(ph1, ph2)
    INTEGER(4), INTENT(IN) :: ph1,ph2
    dllptr(ph1)%chi = ph2
  END SUBROUTINE set_chi_dll

  SUBROUTINE set_nchi_sll(ph1, ph2)
    INTEGER(4), INTENT(IN) :: ph1,ph2
    sllptr(ph1)%nchi = ph2
  END SUBROUTINE set_nchi_sll

  SUBROUTINE set_dchi_sll(ph1, ph2)
    INTEGER(4), INTENT(IN) :: ph1,ph2
    sllptr(ph1)%dchi = ph2
  END SUBROUTINE set_dchi_sll
!----------------------------------------------
  SUBROUTINE deallocate_uf
    deallocate(djsdata,ptrt)
    deallocate(sllptr)
    deallocate(dllptr)
  END SUBROUTINE deallocate_uf
  
   INTEGER(4) FUNCTION find_uf(ph)
    INTEGER(4),INTENT(IN) :: ph
    INTEGER(4) :: parent
    integer(4) :: ptr
    ptr=ph
    if(ptr.eq.djsdata(ph)%par) then
       find_uf=ptr       
       return
    else
       parent=djsdata(ptr)%par
       DO WHILE (ptr /= parent )
         ptr=parent          
         parent=djsdata(parent)%par
      ENDDO
      find_uf=ptr  
       RETURN
    endif
  END FUNCTION find_uf  

END MODULE OBJECT_MODULE3
