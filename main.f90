!--------------------!
module configuration
!----------------------------------------------!
! Most important parameters and data structures
!----------------------------------------------!
save

integer :: lx       ! system length in x direction
integer :: nn       ! number of sites;
integer :: nb       ! number of bonds (depends on lx,ly and boundary conditions)

integer :: mm       ! maximum time length 
integer :: msteps
integer :: nmsr
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

real(8),allocatable :: density(:,:)

!--------------------------!
end module measurementdata
!--------------------------!

!============================!
program dual_unitary_simulator
!=====================================================!
!------------------------------------------------------
use configuration; implicit none

integer :: i,j,nbins,im
!include 'mpif.h'
integer :: ierr, num_process, myrank,root_process
!integer :: status(MPI_STATUS_SIZE),an_id
!call MPI_INIT(ierr)
!call MPI_COMM_RANK (MPI_COMM_WORLD, myrank, ierr)
!call MPI_COMM_SIZE (MPI_COMM_WORLD, num_process, ierr)


myrank = 0
open(10,file='input.in',status='old')
read(10,*)lx,mm,nrep,alpha,nbins
close(10)

msteps= 1

mm=lx/2+1

if (mod(lx,4)/=0) then
    print*,"caution,head on left"
endif

call initran(1,myrank)   
call makelattice()
call initconfig()

   
if (myrank==0) then
   print*,'===================Starting==================='
   print*,'1d pauli string simulator'
   print*,"lx=",lx
   print*,"alpha=",alpha
   print*,"nbins=",nbins
   do im=mm,1,-1
     print*,"init",spin(1,im,:)
   enddo
   print*,'===================Starting==================='
endif

do j=1,nbins
    if (myrank==0) then
        print*,"bin",j
    endif
    
    nmsr=0
   do i=1,msteps
      call update()
      call measure_density()
   enddo
   call writeresults(myrank)
enddo
   
call deallocateall()

end program dual_unitary_simulator
!================================!

!---------------------------!
subroutine update()
!------------------------------------------!
! Carries out one sweep of diagonal updates
!------------------------------------------!
use configuration; implicit none

integer :: i,b,op,s1,s2,irep
real(8), external :: ran

do irep=1,nrep
    do i=1,mm-1
        do b=1+mod(i,2),nb,2
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
enddo

end subroutine update
!-----------------------------!
    
!----------------------------!
subroutine measure_density()
use configuration; use measurementdata; implicit none

integer :: im,s,s1_here,s2_here

do im=1,mm
    do s=1,nn
        density(im,s)=density(im,s)+spin(1,im,s)
    enddo
enddo

nmsr=nmsr+1

end subroutine measure_density

!----------------------!

!-------------------------------!
subroutine writeresults(myrank)
use configuration; use measurementdata; implicit none

integer :: myrank,im
character(len=20)::cha

density = density/dfloat(nmsr)

write(cha,'(i20)')myrank
open(10,file='res'//trim(adjustl(cha))//'.dat',status='replace',position='append')
    write(10,'(4f15.7)') density(:,:)
close(10)

if (myrank==0) then
    do im=mm,1,-1
      print*,"rho",density(im,:)
    enddo
endif

density(:,:)=0.d0
nmsr=0
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


allocate(spin(nrep,mm,nn))

spin(:,:,:)=0
spin(:,1,nn/2)=1


allocate (density(mm,nn))
density(:,:) = 0.d0

end subroutine initconfig
!-------------------------!

!------------------------!
subroutine makelattice()

   use configuration; implicit none
   
   integer :: s,x1,x2,y1,y2,z1,z2
   
   nn=lx
   nb = nn
   allocate(bsites(2,nb))

   do s=1,lx
       bsites(1,s)=s
       bsites(2,s)=s+1
   enddo
   bsites(2,nn)=1
   
   end subroutine makelattice
   !--------------------------!

!--------------------------!
subroutine deallocateall()
!--------------------------! 
use configuration; use measurementdata; implicit none

deallocate (spin)
deallocate (bsites)
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
integer(4) :: w,nb,b
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