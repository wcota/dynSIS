module mod_read_tools
	integer						:: inp_pos
	
contains

	subroutine read_i(i1,a1)
		integer :: i1
		character(len=*) :: a1
		
		write(0,'(a)') 'read: '//trim(adjustl(a1))
		read(*,*) i1;
		write(*,'(a,i)') 'r# '//trim(adjustl(a1))//' = ', i1 
	end subroutine	

	subroutine read_l(l1,a1)
		logical :: l1
		integer :: tmp
		character(len=*) :: a1
		
		write(0,'(a)') 'read: '//trim(adjustl(a1))
		read(*,*) tmp
		l1 = .false.
		if (tmp == 1) l1 = .true.
		write(*,'(a,l)') 'r# '//trim(adjustl(a1))//' = ', l1
	end subroutine
	
	subroutine read_f(f1,a1)
		real*8 :: f1
		character(len=*) :: a1
		
		write(0,'(a)') 'read: '//trim(adjustl(a1))
		read(*,*) f1;
		write(*,'(a,f)') 'r# '//trim(adjustl(a1))//' = ', f1 
	end subroutine
	
	subroutine read_a(a1,a2)
		real*8 :: f1
		character(len=*) :: a1, a2
		
		write(0,'(a)') 'read: '//trim(adjustl(a2))
		read(*,*) a1
		write(*,'(a,a)') 'r# '//trim(adjustl(a2))//' = ', trim(adjustl(a1))
	end subroutine
	
	subroutine read_arg(a1)
		character(len=*) :: a1
		
		call getarg(inp_pos,a1)
		inp_pos = inp_pos + 1
	end subroutine
	
	! http://computer-programming-forum.com/49-fortran/3d6ba18f17fbf6af.htm
	function strip(s1)  result (s2) 
		character(*) :: s1 
		character(char_count(s1)) :: s2 
		integer :: i, n 
		n = 0 
		do i = 1,len_trim(s1)
		   if (s1(i:i) == ' ') cycle 
		   n = n+1 
		   s2(n:n) = s1(i:i) 
		end do 
	end function 
	
	pure function char_count(s)  result (n) 
		character(*),intent(in) :: s 
		integer :: n, i 
		n = 0 
		do i = 1,len_trim(s) 
		   if (s(i:i) == ' ') cycle 
		   n = n+1 
		end do 
	end function 
end module