module UTIL_TEXTURE_MODULE
  implicit none

contains

  subroutine MAKE_GLCM_FEATURES_1D(refdata,adjdata,length,output,nbin,scale)
    implicit none
    integer(4) :: i,j,k,ii,jj
    integer(4),intent(in) :: length
    integer(4),intent(in) :: nbin,scale
    real(4),intent(in) :: refdata(length),adjdata(length)
    real(4),intent(out) :: output
    integer(4) :: refdata2(length),adjdata2(length)
    real(4) :: glcm(0:nbin,0:nbin)
    real(4) :: homogeneity
    refdata2(:)=nint(refdata(:)*scale)
    adjdata2(:)=nint(adjdata(:)*scale)
    glcm(:,:)=0.

    where(refdata2.gt.nbin) refdata2=nbin
    where(adjdata2.gt.nbin) adjdata2=nbin

    do i=1,length
       ii=refdata2(i)
       jj=adjdata2(i)
       glcm(ii,jj)=1.+glcm(ii,jj)
    enddo

    glcm(:,:)=glcm(:,:)/real(length)

    homogeneity=0

    do j =0,nbin
       do i =0,nbin
          homogeneity=glcm(i,j)/real(1+abs(i-j))+homogeneity
       enddo
    enddo

    output=homogeneity
!    print *,output
  endsubroutine MAKE_GLCM_FEATURES_1D

end module Util_Texture_Module
