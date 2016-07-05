! File: dynamics.f90
! ## Running dynamics ##
! ## See README.md to more information and use ##
!-----------------------------------------------------------------------------
! SIS epidemic model algorithm based on the article "Simulation of Markovian epidemic models on large networks"
! Copyright (C) 2016 Wesley Cota, Silvio C. Ferreira
! 
! Please consider citing our article if this code is used.
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
!-----------------------------------------------------------------------------
! Author	: Wesley Cota
! Email		: wesley.cota@ufv.br
! Date		: June 2016
!-----------------------------------------------------------------------------
! see README.md for more details

program mark_dynamics
use mod_netdata
use mod_random
use mod_read_tools
implicit none

	character*1024				:: f_savepath
	
	! Network structure variables
	integer						:: net_kmax
	
	! Dynamics samples variables
	integer						:: dynp_i, dynp_ini, dynp_fin
	real*8						:: dynp_lb
	real*8						:: dynp_save_dt ! uniform for the beginning
	real*8						:: dynp_tmax
	real*8						:: dynp_pINI ! porcentage of network first infected (random)
	logical						:: dynp_new_ini ! redo initial condition randomly for each sample
	
	
	! Dynamics variables
	real*8						:: dyn_p, dyn_q, dyn_dt, dyn_t
	integer, allocatable		:: dyn_ocp(:), dyn_sig(:)
	integer						:: dyn_voc
	integer, allocatable		:: dyn_ini_ocp(:), dyn_ini_sig(:)
	integer						:: dyn_ini_voc, dyn_ini_sk
	integer						:: dyn_sk ! edges from infected vertices
	
	real*8						:: dyn_dt_pos
	
	inp_pos = 1 ! initial value of aux to read command arguments
	
	call read_arg(f_input)
	call read_arg(f_savepath)
	call readEdges()
	net_kmax = maxval(net_k)
	
	call random_ini(iseed)
	
	! We are ready! All Network data is here, now we need the dynamics.
	call read_dyn_parameters()
	
	open(1,file=strip(f_savepath//'/inf_vs_time.dat'))
	
	do dynp_i=dynp_ini,dynp_fin
		call random_select_seed(dynp_i)
		call read_dyn_run()
		write(1,*) ""
	enddo
	
contains

	subroutine read_dyn_run()
	
		integer :: pos_ocp
		integer :: ver, pos_nei
		real*8 :: rnd
	
		if (dynp_new_ini) then
			call random_initial_condition()
		else
			call revert_initial_condition()
		end if
		
		dyn_t = 0d0
		dyn_dt_pos = dynp_save_dt
		dyn_time_loop : do while (dyn_t <= dynp_tmax)
			
			dyn_p = 1d0*dyn_voc/ (dyn_voc + 1d0*dynp_lb * dyn_sk)
			rnd = random_d()
			
			if (rnd < dyn_p) then
				pos_ocp = random_int(1,dyn_voc)
				ver = dyn_ocp(pos_ocp)
				
				! healed
				dyn_sig(ver) = 0
				dyn_sk = dyn_sk - net_k(ver)
				dyn_ocp(pos_ocp) = dyn_ocp(dyn_voc) ! change positions
				dyn_voc = dyn_voc - 1
			else
				select_neig : do
					pos_ocp = random_int(1,dyn_voc)
					ver = dyn_ocp(pos_ocp)
					if (random_d() < 1d0*net_k(ver)/(1d0*net_kmax)) exit select_neig
				enddo select_neig
				
				pos_nei = random_int(net_ini(ver) , net_ini(ver) + net_k(ver) - 1)
				ver = net_con(pos_nei)
				
				if (dyn_sig(ver) == 0) then
					dyn_sig(ver) = 1
					dyn_sk = dyn_sk + net_k(ver)
					dyn_voc = dyn_voc + 1
					dyn_ocp(dyn_voc) = ver
				endif
			endif
			
			if (dyn_voc == 0) exit dyn_time_loop
			
			dyn_dt = 1d0/(dyn_voc + 1d0*dynp_lb * dyn_sk)
			
			dyn_t = dyn_t + dyn_dt
			
			if (dyn_t > dyn_dt_pos) then
				write(1,*) dyn_t, dyn_voc
				dyn_dt_pos = dyn_dt_pos + dynp_save_dt
			endif
			
		enddo dyn_time_loop
		
		
	end subroutine
	
	subroutine random_initial_condition()
		integer :: ver, vti
		
		dyn_sig = 0
		dyn_ocp = 0
		dyn_voc = 0
		dyn_sk = 0
		do vti = 1, net_N*dynp_pINI
			vti_ver : do
				ver = random_int(1,net_N)
				if (dyn_sig(ver) == 0) then
					dyn_voc = dyn_voc + 1
					dyn_ocp(dyn_voc) = ver
					dyn_sig(ver) = 1
					dyn_sk = dyn_sk + net_k(ver)
					exit vti_ver
				endif
			enddo vti_ver
		enddo
		
		if (.not. dynp_new_ini) then
			dyn_ini_ocp = dyn_ocp
			dyn_ini_voc = dyn_voc
			dyn_ini_sig = dyn_sig ! remover depois!
			dyn_ini_sk = dyn_ini_sk
		endif
	end subroutine

	subroutine revert_initial_condition()
		dyn_ocp = dyn_ini_ocp
		dyn_voc = dyn_ini_voc
		dyn_sig = dyn_ini_sig
		dyn_sk = dyn_ini_sk
	end subroutine

	subroutine read_dyn_parameters()
		call read_i(dynp_ini,"Dynamic sample first ID")
		call read_i(dynp_fin,"Dynamic sample last ID")
		call read_f(dynp_lb,"Value of lambda, infection rate (mu is defined as 1)")
		call read_f(dynp_tmax,"Maximum time steps (it stops if the absorbing state is reached)")
		call read_f(dynp_save_dt,"Interval of steps to save data (warning: file can be big!)")
		call read_f(dynp_pINI,"Porcentage of infected vertex on the network.")
		
		dynp_pINI = 1d0*dynp_pINI/100d0
		
		call read_l(dynp_new_ini,"The initial condition must be the different for each sample? 0) to no; 1) to yes.")
		
		
		allocate(dyn_ocp(net_N),dyn_sig(net_N))
		if (.not. dynp_new_ini) then
			allocate(dyn_ini_ocp(net_N),dyn_ini_sig(net_N))
		end if
		
		call random_initial_condition()
	end subroutine
	
end program