module util_shadow_module
  !  use basic_image_analysis_module
  implicit none

contains

  subroutine MAKE_SHADOW_PROB_MAP_1(image,width,length,band,output,alpha,beta,gamma,tau)
    !Automatic and accurate shadow detection using near-infrared information by Rufenacht [2014]
    !original numbers for alpha, beta and gamma are 14, 0.5 and 2,2, respectively.
    implicit none
    integer(4) :: i,j,k,ii,jj
    integer(4),intent(in) :: width,length,band
    real(4),intent(in) :: alpha,beta,gamma,tau
    real(4),intent(in) :: image(width,length,band)
    real(4),intent(out) :: output(width,length)
    real(4) ::brightness(width,length),nlf_v(width,length),nlf_n(width,length)
    real(4) :: candidate(width,length)
    real(4) :: ratio(width,length,3),NIRratio(width,length)
    real(4) :: shadow_map(width,length)
    print *,'calc shadow candidate'
    print *,width,length,band
    brightness(:,:)=(image(:,:,1)+image(:,:,2)+image(:,:,3))/3.
    print *,'brightness'
    nlf_v(:,:)=1/(1+exp(-alpha*(1-brightness(:,:)**(1/gamma)-beta)))
    nlf_n(:,:)=1/(1+exp(-alpha*(1-image(:,:,4)**(1/gamma)-beta)))
    print *,'candidate'
    candidate(:,:)=nlf_v(:,:)*nlf_n(:,:)

    ratio(:,:,1)=image(:,:,1)/image(:,:,4)
    ratio(:,:,2)=image(:,:,2)/image(:,:,4)
    ratio(:,:,3)=image(:,:,3)/image(:,:,4)
    print *,'NIRratio'
    do j =1,length
       do i=1,width
       NIRratio(i,j)=min(maxval(ratio(i,j,:)),tau)/tau
    enddo
 enddo
    print *,'shadomap'

    shadow_map(:,:)=(1-candidate(:,:))*(1-NIRratio(:,:))

    output(:,:)=shadow_map(:,:)
    
  endsubroutine MAKE_SHADOW_PROB_MAP_1

  subroutine MAKE_SHADOW_PROB_MAP_2(image,width,length,band,output,lambda)
    !Correction of shadowing in imaging spectroscopy data by quantification
    !of the proportion of diffuse illumination by schlapfer [2013]    
    implicit none
    integer(4) :: i,j,k,ii,jj
    real(4) :: b
    integer(4),intent(in) :: width,length,band
    real(4),intent(in) :: image(width,length,band)
    real(4),intent(in) :: lambda(band)
    real(4),intent(out) :: output(width,length,3)
    real(4) :: RBindex(width,length)      ! shading effect considering vegetation.
    real(4) :: GRBindex(width,length)     ! shadows on surface types which are not well done by the RBindex
    real(4) :: RNindex(width,length)           ! very dark shadows
    real(4) :: candidate(width,length)
    real(4) :: ratio(width,length,3),NIRratio(width,length)
    real(4) :: shadow_map(width,length)

    do j = 1, length
       do i =1,width
          RNindex(i,j)=min(image(i,j,3),image(i,j,4))
       enddo
    enddo

    RBindex(:,:)=0.7*(image(:,:,3)+0.1*image(:,:,4))/image(:,:,1)

    do j =1,length
       do i=1,width
          b=(image(i,j,3)-image(i,j,2))/(lambda(3)-lambda(2))
          GRBindex(i,j)=1.25*(image(i,j,2)-b*(lambda(2)-lambda(1)))/image(i,j,1)-0.8
       enddo
    enddo

    !step 1
    do j=1,length
       do i=1,width
          if(RNindex(i,j).le.0.025) then
             RBindex(i,j)=RBindex(i,j)*RNindex(i,j)/0.025
          endif
       enddo
    enddo
    
    output(:,:,1)=RNindex(:,:)
    !NIRratio(:,:)
    output(:,:,2)=RBindex(:,:)
    !candidate(:,:)
    output(:,:,3)=GRBindex(:,:)
    
  endsubroutine MAKE_SHADOW_PROB_MAP_2
  
  subroutine MAKE_SHADOW_PROB_MAP_3(image,width,length,band,output,alpha,beta,gamma,tau)
    !New algorithm developed by K. Nakata    
    implicit none
    integer(4) :: i,j,k,ii,jj
    real(4) :: b
    real(4),intent(in) :: alpha,beta,gamma,tau    
    integer(4),intent(in) :: width,length,band
    real(4),intent(in) :: image(width,length,band)
    real(4),intent(out) :: output(width,length,3)
    real(4) ::brightness(width,length),nlf_v(width,length),nlf_n(width,length)
    real(4) :: RBindex(width,length)      ! shading effect considering vegetation.
    real(4) :: candidate(width,length)
    real(4) :: ratio(width,length,3),NIRratio(width,length)
    real(4) :: shadow_map(width,length)

    brightness(:,:)=(image(:,:,1)+image(:,:,2)+image(:,:,3))/3.

    nlf_v(:,:)=1/(1+exp(-alpha*(1-brightness(:,:)**(1/gamma)-beta)))
    nlf_n(:,:)=1/(1+exp(-alpha*(1-image(:,:,4)**(1/gamma)-beta)))
    
    candidate(:,:)=nlf_v(:,:)*nlf_n(:,:)

    ratio(:,:,1)=image(:,:,1)/image(:,:,4)
    ratio(:,:,2)=image(:,:,2)/image(:,:,4)
    ratio(:,:,3)=image(:,:,3)/image(:,:,4)

    do j =1,length
       do i=1,width
       NIRratio(i,j)=min(maxval(ratio(i,j,:)),tau)/tau
    enddo
 enddo

    RBindex(:,:)=0.7*(image(:,:,3)+0.1*image(:,:,4))/image(:,:,1)
    
    output(:,:,1)=candidate(:,:)
    !NIRratio(:,:)
    output(:,:,2)=RBindex(:,:)
    !candidate(:,:)
    output(:,:,3)=NIRratio(:,:)
    
  endsubroutine MAKE_SHADOW_PROB_MAP_3


  
end module Util_Shadow_Module
