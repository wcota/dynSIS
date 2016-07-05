! ## File: dynamics.f90
! ## - main program: running the SIS dynamics, based on SIS II algorithm.
! ## See README.md for more information and use
!-----------------------------------------------------------------------------
! SIS epidemic model algorithm based on the article "Simulation of Markovian epidemic models on large networks"
! Copyright (C) 2016 Wesley Cota, Silvio C. Ferreira
! 
! Please cite the above cited paper as reference to our code.
! 
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------
! Author    : Wesley Cota
! Email     : wesley.cota@ufv.br
! Date      : July 2016
!-----------------------------------------------------------------------------
! See README.md for more details
! This code is available at <https://github.com/wcota/dynSIS>

program dynSIS
use mod_netdata
use mod_random
use mod_read_tools
implicit none
    
    character*1024                :: f_output
    integer, parameter            :: und_output = 10
    
    ! Network structure variables
    integer                       :: net_kmax
    
    ! Dynamics samples variables
    integer                       :: dynp_i, dynp_sam
    real*8                        :: dynp_lb
    real*8                        :: dynp_tmax
    real*8                        :: dynp_pINI ! fraction of network first infected (random)
    
    ! Dynamics variables
    real*8                        :: dyn_p, dyn_q, dyn_dt, dyn_t
    integer, allocatable          :: dyn_ocp(:), dyn_sig(:)
    integer                       :: dyn_voc
    integer                       :: dyn_sk ! edges from infected vertices
    
    real*8                        :: dyn_dt_pos
    
    ! initial value of input counter to read command arguments, used by mod_read_tools
    inp_pos = 1 
    
    ! Files arguments
    call read_arg(f_input)
    call read_arg(f_output)
    
    ! Read network to memory
    call readEdges()
    
    ! To be used in the algorithm.
    net_kmax = maxval(net_k)
    
    ! Initate the random generator
    call random_ini()
    
    ! We are ready! All Network data is here, now we need the dynamics parameters.
    call read_dyn_parameters()
    
    ! Open file and let's run the dynamics.
    open(und_output,file=f_output)
    
    do dynp_i=1,dynp_sam
        call read_dyn_run()
        write(und_output,*) ""
    enddo
    
    close(und_output)
    
    stop 'ALERT: Todo average and write more information.'
    
contains

    subroutine read_dyn_run()
    
        integer :: pos_ocp
        integer :: ver, pos_nei
        real*8  :: rnd
    
        call random_initial_condition()
        
        dyn_t = 0d0
        dyn_dt_pos = 1d0
        
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
                write(und_output,*) dyn_t, 1d0*dyn_voc/net_N
                dyn_dt_pos = dyn_dt_pos + 1d0
            endif
            
        enddo dyn_time_loop
        
    end subroutine
    
    subroutine random_initial_condition()
        integer :: ver, vti
        
        dyn_sig = 0
        dyn_ocp = 0
        dyn_voc = 0
        dyn_sk = 0
        do vti = 1, int(net_N*dynp_pINI)
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
    end subroutine

    subroutine read_dyn_parameters()
        call read_i(dynp_sam,"How much dynamics samples?")
        call read_f(dynp_lb,"Value of infection rate lambda (mu is defined as equal to 1)")
        call read_f(dynp_tmax,"Maximum time steps (it stops if the absorbing state is reached)")
        call read_f(dynp_pINI,"Fraction of infected vertices on the network as initial condition (is random to each sample)")
        
        allocate(dyn_ocp(net_N),dyn_sig(net_N))
    end subroutine
    
end program
