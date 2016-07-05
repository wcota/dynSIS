module mod_random
	
	! Random generator
	integer, parameter :: iseed = 378592914
	! \ Change if you want
	
contains
	
! Random things
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Random number generator KISS05 after a suggestion by George Marsaglia
! in "Random numbers for C: The END?" posted on sci.crypt.random-numbers
! in 1999
!
! version as in "double precision RNGs" in  sci.math.num-analysis  
! http://sci.tech-archive.net/Archive/sci.math.num-analysis/2005-11/msg00352.html
!
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
! Overall period > 2^123  
! 
! 
! A call to rkiss05() gives one random real in the interval [0,1),
! i.e., 0 <= rkiss05 < 1
!
! Before using rkiss05 call kissinit(seed) to initialize
! the generator by random integers produced by Park/Millers
! minimal standard LCG.
! Seed should be any positive integer.
! 
! FORTRAN implementation by Thomas Vojta, vojta@mst.edu
! built on a module found at www.fortran.com
! 
! 
! History:
!        v0.9     Dec 11, 2010    first implementation
!        V0.91    Dec 11, 2010    inlined internal function for the SR component
!        v0.92    Dec 13, 2010    extra shuffle of seed in kissinit 
!        v093     Aug 13, 2012    changed inter representation test to avoid data statements
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	FUNCTION random_d() ! eh o KISS na verdade
		implicit none

		integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
		integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 
		real(r8b),parameter    :: am=4.656612873077392578d-10       ! multiplier 1/2^31

		real(r8b)             :: random_d
		integer(i4b)          :: kiss
		integer(i4b)          :: x,y,z,w              ! working variables for the four generators
		common /kisscom/x,y,z,w 

		x = 69069 * x + 1327217885
		y = ieor (y, ishft (y, 13)); y= ieor (y, ishft (y, -17)); y= ieor (y, ishft (y, 5))
		z = 18000 * iand (z, 65535) + ishft (z, - 16)
		w = 30903 * iand (w, 65535) + ishft (w, - 16)
		kiss = ishft(x + y + ishft (z, 16) + w , -1)
		random_d=kiss*am
	END FUNCTION

	SUBROUTINE random_ini(iinit)
		implicit none
		integer,parameter     :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
		integer,parameter     :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

		integer(i4b) idum,ia,im,iq,ir,iinit
		integer(i4b) k,x,y,z,w,c1,c2,c3,c4
		real(r8b)    rdum
		parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
		common /kisscom/x,y,z,w

		!!! Test integer representation !!!
		c1=-8
		c1=ishftc(c1,-3)
		!     print *,c1
		if (c1.ne.536870911) then
			print *,'Nonstandard integer representation. Stoped.'
			stop
		endif

		idum=iinit
		idum= abs(1099087573 * idum)               ! 32-bit LCG to shuffle seeds
		if (idum.eq.0) idum=1
		if (idum.ge.IM) idum=IM-1

		k=(idum)/IQ
		idum=IA*(idum-k*IQ)-IR*k
		if (idum.lt.0) idum = idum + IM
		if (idum.lt.1) then
				x=idum+1 
			else 
				x=idum
		endif
		k=(idum)/IQ
		idum=IA*(idum-k*IQ)-IR*k
		if (idum.lt.0) idum = idum + IM
		if (idum.lt.1) then 
				y=idum+1 
			else 
				y=idum
		endif
		k=(idum)/IQ
		idum=IA*(idum-k*IQ)-IR*k
		if (idum.lt.0) idum = idum + IM
		if (idum.lt.1) then
				z=idum+1 
		else 
				z=idum
		endif
		k=(idum)/IQ
		idum=IA*(idum-k*IQ)-IR*k
		if (idum.lt.0) idum = idum + IM
		if (idum.lt.1) then
				w=idum+1 
			else 
		w=idum
		endif

		rdum=random_d()

		return
      end subroutine
	
	function random_int(i1,i2)
		integer, intent(in) :: i1, i2
		integer :: random_int
		
		random_int = min(int(random_d()*(i2+1-i1))+i1,i2)
	end function
	
	subroutine random_select_seed(i1)
		integer, intent(in) :: i1
		integer :: iran
		real*8 :: ran
		
		call random_ini(iseed)
		do iran=1,100000*i1
			ran = random_d()
		enddo	
	end subroutine
	
	! / Random things
	
end module