! File: mod_netdata.f90
! ## Conversion of list of edges from ascii to binary ##
! ## See README.md to more information and use ##
!-----------------------------------------------------------------------------
! SIS epidemic model algorithm based on the article "Simulation of Markovian epidemic models on large networks"
! Copyright (C) 2016 Wesley Cota, Silvio C. Ferreira
! 
! This code is a free software: you can redistribute it and/or modify it under the terms 
! of the GNU Lesser General Public License as published by the Free Software Foundation, 
! either version 3 of the License, or (at your option) any later version.
! 
! It is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
! 
! GNU Lesser General Public License is available at <http://www.gnu.org/licenses/>.
!
! Use: none. This is just a module to be used in another program.
!-----------------------------------------------------------------------------
! Author	: Wesley Cota
! Email		: wesley.cota@ufv.br
! Date		: June 2016
!-----------------------------------------------------------------------------
! see README.md for more details

module mod_netdata
implicit none

	! file input names and signals
	character*1024				:: f_input
	integer						:: f_sig ! file read signal
	
	! vertices and edges
	integer						:: net_N ! total vertices of the network
	integer						:: net_skk ! total number of edges
	
	! network data and adjacency list
	integer, allocatable		:: net_con(:), net_ini(:), net_k(:) ! adjacency, edges positions, auxiliar edges pos., degree of each vertex lists
	
contains
	subroutine readEdges()
	
		! input file read variables
		integer						:: f_input_nl ! number of lines
		
		integer, allocatable		:: tmp_con(:,:) ! temporary edges lists
		integer						:: pos_con ! position of edges and vertices
		integer 					:: ver, ver1, ver2
		
		integer, allocatable		:: aux_ini(:)
		
		! auxiliar integers
		integer						:: aux_i1, aux_i2
		integer						:: iaux, itmp
		integer						:: i, j
			
		! input file is opened
		open(1,file=f_input,action='read',status='old')
		
		! We set net_N to be the largest value read
		net_N = 0
	
		! calculating number of lines to read and allocate matrices...
		f_input_nl = 0
		f_input_loop : do
			read(1,*,IOSTAT=f_sig) aux_i1, aux_i2
			if (f_sig < 0) exit f_input_loop ! if end of file, exit. No "goto" allowed! :)
			
			if (aux_i1 == aux_i2) call print_error("Self-connection found! Verify your data.") ! Self connection. Something is wrong with data!
			if (aux_i1 < 1 .or. aux_i2 < 1) call print_error('Vertex id MUST be >= 1. Verify your data!')
			
			net_N = max(net_N, aux_i1, aux_i2) ! will be the largest value found
			
			f_input_nl = f_input_nl + 1
		enddo f_input_loop
	
		close(1); open(1,file=f_input,action='read',status='old') ! input file is opened again

		! We already have the size of the network net_N. Now, assuming all edges undirected, net_skk = 2*f_input_nl, twice the number of lines (see README.md to know how to save the data).
		net_skk = 2*f_input_nl
	
		allocate(tmp_con(2,f_input_nl)) ! We need to calculate the degree of the network, and to do this we will have to read the file again. So, we save the data on the way.
		allocate(net_k(net_N))
	
		! number of edges added
		pos_con = 0
		tmp_con = 0 ! list of connections read
		net_k = 0 ! degree is zero at the beginning
	
		! Let's read the file again, this time saving the connections, since we already know the total number of connections.
		input_con_loop : do
			read(1,*,IOSTAT=f_sig) aux_i1, aux_i2
			if (f_sig < 0) exit input_con_loop ! if end of file, exit.

			pos_con = pos_con + 1 ! new connection is saved
		
			! save line to array
			tmp_con(1,pos_con) = aux_i1
			tmp_con(2,pos_con) = aux_i2
		
			! degree is added to both vertices
			net_k(aux_i1) = net_k(aux_i1) + 1
			net_k(aux_i2) = net_k(aux_i2) + 1
		enddo input_con_loop
		close(1) ! we do not need the file anymore
	
		! Is REALLY everything ok? Let's check!
		if (count(tmp_con == 0) > 0) call print_error('Please, verify your data. Something is not right! [count(tmp_con == 0) > 0]')
		if (sum(net_k) .ne. net_skk) call print_error('Please, CHECK! Sum of degrees IS NOT equal the number of connections found & 
in the file')
	
		! Everything OK! 
		! Now we will build the REAL list of adjacency.
		call make_ini() ! This will use the degrees array and build the net_ini matrix.
		allocate(aux_ini(net_N)) ! We will use to add the valid connections. It tell us the last *free* position to add to the net_con list.
	
		! To begin, we set the initial values.
		aux_ini = net_ini
	
		! We allocate the real list of adjacency.
		allocate(net_con(net_skk))
	
		! Now we read again the data and add the connections
		net_con = 0 ! We do not have any connection added.
		do iaux=1,f_input_nl
			ver1 = tmp_con(1,iaux)
			ver2 = tmp_con(2,iaux)
			
			! Just to check...
			if (net_con(aux_ini(ver1)) .ne. 0) call print_error('Houston, we have had a problem. Verify your data! &
[net_con(aux_ini(ver1)) .ne. 0]')
			if (net_con(aux_ini(ver2)) .ne. 0) call print_error('Houston, we have had a problem. Verify your data! &
[net_con(aux_ini(ver2)) .ne. 0]')
			
			! Ok...
			net_con(aux_ini(ver1)) = ver2
			aux_ini(ver1) = aux_ini(ver1) + 1
			
			net_con(aux_ini(ver2)) = ver1
			aux_ini(ver2) = aux_ini(ver2) + 1
		enddo
		
		! Let's check again...
		if (count(net_con == 0) > 0) call print_error('Verify your data. [count(net_con == 0) > 0]')
	
		! We do not need tmp_con anymore!
		deallocate(tmp_con)
		
		! Now, everything is fine and we have the network data ready to be used! ;)
	
	contains

		subroutine make_ini() 
			
			call deal(net_ini)
			allocate(net_ini(net_N))
			net_ini = 0
	
			iaux = 1
			do ver=1,net_N
				itmp = net_k(ver)
				net_ini(ver) = iaux
				iaux = iaux + itmp
			enddo
		end subroutine 
	
	end subroutine readEdges

	
	! Common and useful subroutines
	subroutine deal(arr)
		integer,allocatable :: arr(:)
		
		if (allocated(arr)) deallocate(arr)
	end subroutine

	subroutine print_warning(c1)
		character(*) :: c1
		
		write(0,'(a,a)') "##! Alert !## ", c1
	end subroutine

	subroutine print_error(c1)
		character(*) :: c1
		
		write(0,'(a,a)') "##!!! ERROR !!!## ", c1
		stop ""
	end subroutine
	
end module