!--------------------!
module configuration
!----------------------------------------------!
! Most important parameters and data structures
!----------------------------------------------!
save

integer :: lx       ! system length in x direction
integer :: tstep    ! steps between ajacent measurement

integer :: mm       ! maximum time length 
integer :: msteps
integer :: nrep     ! number of replica, renyi index

real(8) :: alpha

integer, allocatable :: spin(:,:,:)      ! spin state
integer, allocatable :: bsites(:,:)  ! list of sites bsites(1,b),bsites(2,b) at bond b

end module configuration
!------------------------!

!----------------------!
module measurementdata
!----------------------------------------------!
! Data were measured quantities are accumulated
! - See 'measure' for explanation of variables
!----------------------------------------------!
save

real(8), allocatable :: purity(:,:)
real(8),allocatable :: density(:,:,:)


!--------------------------!
end module measurementdata
!--------------------------!

!============================!
program dual_unitary_simulator
!=====================================================!
!------------------------------------------------------
use configuration; implicit none

integer :: i,j,nbins,im
integer :: ierr, num_process, myrank,root_process
!include 'mpif.h'
!integer :: status(MPI_STATUS_SIZE),an_id
!call MPI_INIT(ierr)
!call MPI_COMM_RANK (MPI_COMM_WORLD, myrank, ierr)
!call MPI_COMM_SIZE (MPI_COMM_WORLD, num_process, ierr)


myrank = 0
open(10,file='input.in',status='old')
read(10,*)lx,nrep,alpha,nbins
close(10)

msteps=50000

do alpha=0.1,1.01,0.1
    
mm=lx/2+1
tstep = 10

if (mod(lx,4)/=0) then
    print*,"caution,head on left"
    stop
endif

call initran(1,myrank)   
call makelattice()
call initconfig()
   
if (myrank==0) then
   print*,'===================Starting==================='
   print*,'1d pauli string simulator,ResumEE'
   print*,"lx=",lx
   print*,"alpha=",alpha
   print*,"mm=",mm
   print*,"nbins=",nbins
   !do im=mm,1,-1
   !  print*,"init",spin(1,im,:)
   !enddo
   print*,"lx",lx
   print*,'===================Starting==================='
endif

do j=1,nbins
    if (myrank==0) then
        print*,"bin",j
    endif
   call update(0)
   do i=1,msteps
      call update(1)
      call measure_density()
   enddo
   call writeresults(myrank)
enddo
   
call deallocateall()

enddo

end program dual_unitary_simulator
!================================!

!---------------------------!
subroutine update(ismeasure)
!------------------------------------------!
! Carries out one sweep of diagonal updates
!------------------------------------------!
use configuration; implicit none

integer :: i,b,op,s1,s2,irep,tm,ismeasure
real(8), external :: ran

do irep=1,nrep
    do i=1,mm-2 ! update till the second last time slice
        do b=mm-i,lx-(mm-i),2 !only consider bond in side light-cone
            s1=bsites(1,b)
            s2=bsites(2,b)
            if (spin(irep,i,s1)==0.and.spin(irep,i,s2)==0) then
                spin(irep,i+1,s1)=0
                spin(irep,i+1,s2)=0
            elseif (spin(irep,i,s1)-spin(irep,i,s2)/=0) then
                if (ran() <= 1.d0-alpha) then
                    spin(irep,i+1,s1)=spin(irep,i,s2)
                    spin(irep,i+1,s2)=spin(irep,i,s1)
                else
                    spin(irep,i+1,s1)=1
                    spin(irep,i+1,s2)=1
                endif
            else
                if (ran() <= alpha*2.d0/3.d0) then
                    if (ran()<=0.5) then
                        spin(irep,i+1,s1)=1
                        spin(irep,i+1,s2)=0
                    else
                        spin(irep,i+1,s1)=0
                        spin(irep,i+1,s2)=1
                    endif
                else
                    spin(irep,i+1,s1)=1
                    spin(irep,i+1,s2)=1
                endif
            endif
        enddo
    enddo
    !  !---do measurement---!
    if (ismeasure==1) then     
        do tm=tstep+1,mm,tstep
            call measure_EE(tm)
        enddo
    endif
enddo


end subroutine update
!-----------------------------!
    
!----------------------------!
subroutine measure_EE(tm)
use configuration; use measurementdata; implicit none

integer :: im,s,s1,s2,irep,b,spin1,spin2,idx,tm
real(8) :: term,enrg
real(8),allocatable :: out_prob(:,:,:)  ! even site case, no single site left, four cases for a gate 

! (0,0)->1; (0,1)->2; (1,0)->3; (1,1)->4

allocate(out_prob(nrep,tm,4))
out_prob(:,:,:)=0.d0


do b=1,tm-1 ! b is the gate index  
    do irep=1,nrep
        s1 = bsites(1,2*b-1+mm-tm)
        s2 = bsites(2,2*b-1+mm-tm)
        spin1=spin(irep,tm-1,s1)
        spin2=spin(irep,tm-1,s2)
        if (spin1==1.and.spin2==0) then !(1,0)->
            out_prob(irep,b,4)=alpha  !(1,1)
            out_prob(irep,b,2)=1.d0 - alpha
        elseif (spin1==0.and.spin2==1) then !(0,1)->
            out_prob(irep,b,4)=alpha
            out_prob(irep,b,3)=1.d0-alpha
        elseif (spin1==1.and.spin2==1) then !(1,1)->
            out_prob(irep,b,2)=alpha/3.d0 !(0,1)
            out_prob(irep,b,3)=alpha/3.d0 !(1,0)
            out_prob(irep,b,4)=1-2.d0*alpha/3.d0 !(1,1)
        else !(0,0)->
            out_prob(irep,b,1)=1.d0
        endif
    enddo
enddo


enrg = 1.d0
idx=0
!do b=tm-1,1,-1 !measure sequence; from head to tail or inverse
do b=1,tm-1,1
    idx=idx+1
    term=0
    term = term+out_prob(1,b,1)*out_prob(2,b,1)
    term = term+out_prob(1,b,2)*out_prob(2,b,2)/3.d0
    term = term+out_prob(1,b,3)*out_prob(2,b,3)/3.d0
    term = term+out_prob(1,b,4)*out_prob(2,b,4)/9.d0
    enrg = enrg*term
    if (term==0) then
        exit
    endif
    purity(tm/tstep,idx)=purity(tm/tstep,idx)+enrg    
enddo

deallocate(out_prob)

end subroutine measure_EE
!----------------------!
!----------------------------!
subroutine measure_density()
use configuration; use measurementdata; implicit none

integer :: im,s,s1_here,s2_here,irep

do im=1,mm
    do s=1,lx
        do irep=1,nrep
            density(irep,im,s)=density(irep,im,s)+spin(irep,im,s)
        enddo
    enddo
enddo


end subroutine measure_density

!----------------------!

!-------------------------------!
subroutine writeresults(myrank)
use configuration; use measurementdata; implicit none

integer :: myrank,im,s,ti,irep
character(len=20)::cha,cha_t,cha_irep

purity(:,:) = purity(:,:)/dfloat(msteps)/dfloat(nrep)

do ti=1,mm/tstep
    write(cha,'(i20)')myrank
    write(cha_t,'(i20)')ti*tstep
    open(10,file="t"//trim(adjustl(cha_t))//'res'//trim(adjustl(cha))//'.dat',status='unknown',position='append')
    do s=1,ti*tstep
        write(10,'(4f19.15)') purity(ti,s)
    enddo
enddo
close(10)

if (myrank==0) then
    do ti=1,mm/tstep
        do s=1,ti*tstep
          print*,"EE at",s,":",-log(purity(ti,s))
        enddo
        print*,"========================================================="
    enddo
endif
purity(:,:) = 0.d0

!---density---!
density = density/dfloat(msteps)
do irep=1,nrep
    write(cha,'(i20)')myrank
    write(cha_irep,'(i20)')irep
    open(10,file="rep"//trim(adjustl(cha_irep))//'rho'//trim(adjustl(cha))//'.dat',status='replace',position='append')
    write(10,'(4f15.7)') density(irep,:,:)
enddo
close(10)
density(:,:,:)=0.d0

end subroutine writeresults
!---------------------------!

!-----------------------!
subroutine initconfig()
!--------------------------------------------------------------!
! Allocates the most important arrays an initializes the stored
! spin configuration and the empty operator string.
!--------------------------------------------------------------!
use configuration; use measurementdata; implicit none

integer :: i


allocate(spin(nrep,mm,lx))

spin(:,:,:)=0
spin(:,1,lx/2)=1
!spin(:,1,lx/2+1)=1


allocate(purity(lx/tstep,lx/2)) !save every 20 time, each time save every two sites
purity(:,:) = 0.d0

allocate (density(nrep,mm,lx))
density(:,:,:) = 0.d0

end subroutine initconfig
!-------------------------!

!------------------------!
subroutine makelattice()

   use configuration; implicit none
   
   integer :: s,x1,x2,y1,y2,z1,z2
   
   allocate(bsites(2,lx))

   do s=1,lx
       bsites(1,s)=s
       bsites(2,s)=s+1
   enddo
   bsites(2,lx)=1
   
end subroutine makelattice
!--------------------------!

!--------------------------!
subroutine deallocateall()
!--------------------------! 
use configuration; use measurementdata; implicit none

deallocate (spin)
deallocate (bsites)
deallocate (purity)
deallocate (density)
end subroutine deallocateall
!----------------------------!

!----------------------!
real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
implicit none

real(8)    :: dmu64
integer(8) :: ran64,mul64,add64
common/bran64/dmu64,ran64,mul64,add64

ran64=ran64*mul64+add64
ran=0.5d0+dmu64*dble(ran64)

end function ran
!----------------!

!---------------------!
subroutine initran(w,rank)
!---------------------!
implicit none


integer(8) :: irmax
integer(4) :: w,lx,b
integer :: rank, system_time

real(8)    :: dmu64,rdm
integer(8) :: ran64,mul64,add64,ii,i
common/bran64/dmu64,ran64,mul64,add64


irmax=2_8**31
irmax=2*(irmax**2-1)+1
mul64=2862933555777941757_8
add64=1013904243
dmu64=0.5d0/dble(irmax)

call system_clock(system_time)
ran64=abs(system_time-(rank*1989+2010)*20130928)

end subroutine initran
!----------------------!