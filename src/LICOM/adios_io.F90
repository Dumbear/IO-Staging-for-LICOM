!----------------------------------------------------
!   READ AND WRITE THE BP FILE WITH ADIOS 1.3.1
!   Written by ZOU yinlong and XUE wei 2012.6
!----------------------------------------------------

!   Start the adios
subroutine start_adios
#include <def-undef.h>
use param_mod, only:mytid
use msg_mod, only: mpi_comm_ocn
integer :: adios_err

    call adios_init ("licom2.xml", mpi_comm_ocn, adios_err)
    if(adios_err .ne. 0) write(6,*) 'proc:',mytid, &
                                    'error in adios_init',adios_err
end subroutine start_adios

!   Shutdown the adios
subroutine end_adios
#include <def-undef.h>
use param_mod, only:mytid
integer :: adios_err

    call adios_finalize (mytid, adios_err)
    if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_finalize'

end subroutine end_adios

!   Read the restart file with bp format
!   only use the basic method
!   overlapping area between processes is not stored in file
subroutine inirun_adios_read_noneoverlap
#include <def-undef.h>
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

#ifdef COUP
use shr_msg_mod
use shr_mpi_mod
use shr_sys_mod
use buf_mod, only: t_cpl,s_cpl,q,u_cpl,v_cpl,dhdx,dhdy
use controlocn_mod
#endif

implicit none

#ifdef COUP
!----- local  ------
integer            :: fid    ! nc domain file ID
integer            :: dimid  ! nc dimension id
integer            :: vid    ! nc variable ID
integer            :: rcode  ! nc return code
integer            :: ntim   ! temporary
#endif

character (len=21) :: fname
integer   :: gcnt,var_count,attr_count,tmp,adios_err
integer*8 :: adios_groupsize,adios_totalsize,adios_handle, &
             gh,adios_buf_size,read_bytes,total_bytes
integer*8,dimension(0:2) :: starts_2d, starts_3d, &
                            count_2d, count_3d

    if (mytid==0) then
       write(6,*) "beginning of inirun_adios_read_noneoverlap !"
        open (17,file='rpointer.ocn.adios',form='formatted')
        read (17,'(a21)') fname
        close(17)
    endif 

    call mpi_bcast(fname,21,mpi_character,0,mpi_comm_ocn,adios_err)
 
    total_bytes =0

    call adios_fopen (adios_handle,fname,mpi_comm_self,gcnt,adios_err)
    if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_fopen'

    call adios_gopen (adios_handle, gh,"ssaveins",var_count,attr_count,adios_err)
    if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_gopen'

    call adios_read_var (gh,"/var/number_month",0,1,number_month,read_bytes)
    month=number_month
    total_bytes=total_bytes+read_bytes
    call adios_read_var (gh,"/var/number_day",0,1,number_day,read_bytes )
    total_bytes=total_bytes+read_bytes

    count_2d(0)=imt
    count_2d(1)=jmt

    count_3d(0)=imt
    count_3d(1)=jmt
    count_3d(2)=km

    if(iy == (ny_proc-1)) then
        count_3d(1)= jmt_global-j_start(mytid+1)+1
        count_2d(1)= jmt_global-j_start(mytid+1)+1
    endif

    starts_2d(0)=i_start(mytid+1)-1
    starts_2d(1)=j_start(mytid+1)-jst_global

    starts_3d(0)=i_start(mytid+1)-1
    starts_3d(1)=j_start(mytid+1)-jst_global
    starts_3d(2)=0

    !write(6,*) 'READ proc:',mytid,'I: ',count_2d(0),starts_2d(0)
    !write(6,*) 'READ proc:',mytid,'J: ',count_2d(1),starts_2d(1)

    call adios_read_var(gh,"/var/h0" ,starts_2d,count_2d, &
                        h0(1:count_2d(0),1:count_2d(1)),read_bytes)
    total_bytes=total_bytes+read_bytes
    call adios_read_var(gh,"/var/u"  ,starts_3d,count_3d, &
                        u(1:count_3d(0),1:count_3d(1),1:count_3d(2)),read_bytes)
    total_bytes=total_bytes+read_bytes
    call adios_read_var(gh,"/var/v"  ,starts_3d,count_3d, &
                        v(1:count_3d(0),1:count_3d(1),1:count_3d(2)),read_bytes)
    total_bytes=total_bytes+read_bytes
    call adios_read_var(gh,"/var/at1",starts_3d,count_3d, &
                        at(1:count_3d(0),1:count_3d(1),1:count_3d(2),1),read_bytes)
    total_bytes=total_bytes+read_bytes
    call adios_read_var(gh,"/var/at2",starts_3d,count_3d, &
                        at(1:count_3d(0),1:count_3d(1),1:count_3d(2),2),read_bytes)
    total_bytes=total_bytes+read_bytes
    !write(6,*) 'proc:',mytid,'read total_bytes is',total_bytes

#ifdef COUP
     if (nstart==2) then
            month=(cdate/10000-1)*12+mod(cdate,10000)/100
!M
!$OMP PARALLEL DO PRIVATE (J,I)
            do j=1,jmt
            do i=1,imt
               ! lihuimin, TODO, to be considered
               ! lihuimin, 2012.7.23, coordinate with flux_cpl, ft. yu
               !t_cpl (i,j)  = 273.15+at(i,jmt-j+1,1,1)
               !s_cpl (i,j)  = at(i,jmt-j+1,1,2)*1000.+35.
               t_cpl (i,j)  = 273.15+at(i,j,1,1)
               s_cpl (i,j)  = at(i,j,1,2)*1000.+35.
               ! modi end
               q     (i,j)  = 0.0
               u_cpl (i,j)  = 0.0
               v_cpl (i,j)  = 0.0
               dhdx  (i,j)  = 0.0
               dhdy  (i,j)  = 0.0
            end do
            end do
        else
        call adios_read_var(gh,"/var/t_cpl" ,starts_2d,count_2d, &
                        t_cpl(1:count_2d(0),1:count_2d(1)),read_bytes)
          total_bytes=total_bytes+read_bytes
         call adios_read_var(gh,"/var/s_cpl" ,starts_2d,count_2d, &
                        s_cpl(1:count_2d(0),1:count_2d(1)),read_bytes)
          total_bytes=total_bytes+read_bytes
         call adios_read_var(gh,"/var/u_cpl" ,starts_2d,count_2d, &
                        u_cpl(1:count_2d(0),1:count_2d(1)),read_bytes)
          total_bytes=total_bytes+read_bytes
          call adios_read_var(gh,"/var/v_cpl" ,starts_2d,count_2d, &
                        v_cpl(1:count_2d(0),1:count_2d(1)),read_bytes)
          total_bytes=total_bytes+read_bytes
          call adios_read_var(gh,"/var/dhdx" ,starts_2d,count_2d, &
                        dhdx(1:count_2d(0),1:count_2d(1)),read_bytes)
          total_bytes=total_bytes+read_bytes
          call adios_read_var(gh,"/var/dhdy" ,starts_2d,count_2d, &
                        dhdy(1:count_2d(0),1:count_2d(1)),read_bytes)
          total_bytes=total_bytes+read_bytes
          call adios_read_var(gh,"/var/q" ,starts_2d,count_2d, &
                        q(1:count_2d(0),1:count_2d(1)),read_bytes)
          total_bytes=total_bytes+read_bytes
        endif
#endif
   write(*,*) 'number_month =',number_month,'mon0=',mon0,&
                       'number_day=',number_day,'iday=',iday
    
    call adios_gclose(gh, adios_err)
    if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_gclose'
    call adios_fclose(adios_handle, adios_err)
    if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_fclose'

    if((iy==ny_proc-1).and.((j_start(mytid+1)+jmt-1)>jmt_global)) then
        tmp=jmt_global-j_start(mytid+1)+1+1
        h0(1:imt,tmp:jmt) = 0.0
        u (1:imt,tmp:jmt,1:km) = 0.0
        v (1:imt,tmp:jmt,1:km) = 0.0
        at(1:imt,tmp:jmt,1:km,1:2) = 0.0
    endif

    return

end subroutine inirun_adios_read_noneoverlap

!   Write the output file with bp format
!   overlapping area between processes is not stored to file
subroutine ssaveins_adios_noneoverlap

#include <def-undef.h>
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
#ifdef COUP
!LPF 20120815
use buf_mod, only:t_cpl,s_cpl,u_cpl,v_cpl,dhdx,dhdy,q
!LPF 20120815
#endif
use msg_mod, only: mpi_comm_ocn
use output_mod, only: spval

implicit none

character (len=21) :: fname
integer:: ni_offset,nj_offset,nk_offset, &
          ni_local,nj_local,nk_local, &
          ni_global,nj_global,nk_global, &
          ni_start,nj_start,nk_start, &
          ni_end,nj_end,nk_end, &
          nwmf,adios_err,error_sum
integer*8:: adios_groupsize,adios_totalsize,adios_handle,group_comm
real(r8),dimension(imt,jmt) :: buffer_2d
real(r8),dimension(imt,jmt,km) :: buffer_3d

   ! number_day = iday
    number_day = iday+1
    number_month = month
    if ( number_day  > imd ) then
        number_day = 1
        number_month = month + 1
    end if
    nwmf= iyfm
 
    if(mytid==0) write(*,*)'in instant','iday=',iday,'imd=',imd,'rest_freq=',rest_freq
    
    if (mod(iday,rest_freq) == 0 .or. iday == 1) then !only for low monthly output

        fname(1:8)='fort.22.'
        fname(13:13)='-'
        fname(16:16)='-'
       ! write(fname(14:15),'(i2.2)')mon0
        write(fname(14:15),'(i2.2)') number_month   
        write(fname(9:12),'(i4.4)')nwmf
        write(fname(17:18),'(i2.2)') number_day
        fname(19:21)='.bp'

       if(mytid==0) write(6,*)'in instant',iday,imd,rest_freq

        if(mytid==0) write(6,*) 'begin in ssavecdf_adios_noneoverlap: ',fname
        !write(6,*) 'proc:',mytid,'filename of adios is',fname

 !LPF20140119
          call exch_boundary(h0,1)
          call exch_boundary(u,km)
          call exch_boundary(v,km)
          call exch_boundary(at(1,1,1,1),km)
          call exch_boundary(at(1,1,1,2),km)
          call exch_boundary(ws(1,1,1),1)
          call exch_boundary(su,1)
          call exch_boundary(sv,1)
          call exch_boundary(swv,1)
          call exch_boundary(lwv,1)
          call exch_boundary(sshf,1)
          call exch_boundary(lthf,1)
          call exch_boundary(fresh,1)
 !LPF20140119

       call adios_open (adios_handle,"ssaveins",fname,"w",mpi_comm_ocn,adios_err)
        if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_open'

        ni_global= imt_global
        nj_global= jmt_global-jst_global+1
        nk_global= km
       
       ! if(mytid==0) write(6,*) '!!!!!!!!!!!!!!!!!!!!! adios: imt_global,jmt_global,jst_global,nj_global,km:',imt_global,jmt_global,jst_global,nj_global,km
       ! write(6,*)'!!!!!!!!!!!!!!!!!!! mytid,imt,jmt,num_overlap: ',mytid,imt,jmt,num_overlap
        if( nx_proc == 1 ) then
            ni_local = ni_global
            ni_offset= 0
            ni_start = 1
        else if( ix == 0 ) then
            ni_local = imt - num_overlap/2
            ni_offset= 0
            ni_start = 1
        else if( ix == (nx_proc-1) ) then
            ni_local = imt - num_overlap/2
            ni_offset= i_start(mytid+1)
            ni_start = num_overlap/2+1
        else
            ni_local = imt - num_overlap
            ni_offset= i_start(mytid+1)
            ni_start = num_overlap/2+1
        endif

        if( ny_proc == 1 ) then
            nj_local = nj_global
            nj_offset= 0
            nj_start = 1
        else if( iy == 0 ) then
            nj_local = jmt - num_overlap/2
            nj_offset= 0
            nj_start = 1
        else if( iy == (ny_proc-1) ) then
            nj_local = jmt_global - j_start(mytid+1)
            nj_offset= j_start(mytid+1)-jst_global+1
            nj_start = num_overlap/2+1
        else
            nj_local = jmt - num_overlap
            nj_offset= j_start(mytid+1)-jst_global+1
            nj_start = num_overlap/2+1
        endif

        nk_local = km
        nk_offset= 0
        nk_start = 1

        ni_end = ni_start+ni_local-1
        nj_end = nj_start+nj_local-1
        nk_end = nk_start+nk_local-1

        !write(6,*) 'proc:',mytid,'I: ',ni_local,ni_offset,ni_start
        !write(6,*) 'proc:',mytid,'J: ',nj_local,nj_offset,nj_start
        !write(6,*) 'proc:',mytid,'K: ',nk_local,nk_offset,nk_start

        adios_groupsize= 4+4 + &
                         4*9 + &
                         8*5*(ni_local)*(nj_local)*(nk_local) + &
                         8*8*(ni_local)*(nj_local)
#ifdef COUP
        adios_groupsize= adios_groupsize+8*7*(ni_local)*(nj_local)
#endif
        !write(6,*) 'proc:',mytid,'groupsize of adios is',adios_groupsize,'bytes'

        call adios_group_size (adios_handle,adios_groupsize, &
                               adios_totalsize,adios_err)
        if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_group_size'

        error_sum=0
        call adios_write(adios_handle,"ni_global",ni_global,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nj_global",nj_global,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk_global",nk_global,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"ni_offset",ni_offset,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nj_offset",nj_offset,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk_offset",nk_offset,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"ni_local",ni_local,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nj_local",nj_local,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk_local",nk_local,adios_err)
        error_sum = error_sum+adios_err
        if(error_sum .ne. 0) write(6,*) 'proc:',mytid,'error in writing dimensions'

        error_sum=0
        call adios_write (adios_handle,"number_month",number_month,adios_err)
        error_sum = error_sum+adios_err
        call adios_write (adios_handle,"number_day",number_day,adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=h0(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"h0",buffer_2d(ni_start:ni_end, &
                                               nj_start:nj_end),&
                                               adios_err)
        error_sum = error_sum+adios_err

        buffer_3d(:,:,:)=u(:,:,:)
        if(iy==0) then
        buffer_3d(:,1,:)=spval
        endif
        if(ix==0) then
        buffer_3d(1,:,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_3d(imt,:,:)=spval
        endif
        call adios_write (adios_handle,"u",buffer_3d(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk_end),&
                                             adios_err)
        error_sum = error_sum+adios_err

        buffer_3d(:,:,:)=v(:,:,:)
        if(iy==0) then
        buffer_3d(:,1,:)=spval
        endif
        if(ix==0) then
        buffer_3d(1,:,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_3d(imt,:,:)=spval
        endif
        call adios_write (adios_handle,"v",buffer_3d(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk_end),&
                                             adios_err)
        error_sum = error_sum+adios_err

        buffer_3d(:,:,:)=at(:,:,:,1)
        if(iy==0) then
        buffer_3d(:,1,:)=spval
        endif
        if(ix==0) then
        buffer_3d(1,:,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_3d(imt,:,:)=spval
        endif
        call adios_write (adios_handle,"at1",buffer_3d (ni_start:ni_end,&
                                                nj_start:nj_end, &
                                                nk_start:nk_end), &
                                                adios_err)
        error_sum = error_sum+adios_err

        buffer_3d(:,:,:)=at(:,:,:,2)
        if(iy==0) then
        buffer_3d(:,1,:)=spval
        endif
        if(ix==0) then
        buffer_3d(1,:,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_3d(imt,:,:)=spval
        endif
        call adios_write (adios_handle,"at2",buffer_3d(ni_start:ni_end, &
                                                nj_start:nj_end, &
                                                nk_start:nk_end), &
                                                adios_err)
        error_sum = error_sum+adios_err

        buffer_3d(:,:,1:km)= ws(:,:,1:km)
        if(iy==0) then
        buffer_3d(:,1,:)=spval
        endif
        if(ix==0) then
        buffer_3d(1,:,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_3d(imt,:,:)=spval
        endif
        call adios_write (adios_handle,"ws",buffer_3d(ni_start:ni_end, &
                                               nj_start:nj_end, &
                                               nk_start:nk_end),&
                                               adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=su(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif        
        call adios_write (adios_handle,"su",buffer_2d(ni_start:ni_end, &
                                               nj_start:nj_end),&
                                               adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=sv(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"sv",buffer_2d(ni_start:ni_end, &
                                               nj_start:nj_end),&
                                               adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=swv(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"swv",buffer_2d(ni_start:ni_end, &
                                                 nj_start:nj_end),&
                                                 adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=lwv(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"lwv",buffer_2d(ni_start:ni_end, &
                                                 nj_start:nj_end),&
                                                 adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=sshf(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"sshf",buffer_2d(ni_start:ni_end, &
                                                   nj_start:nj_end),&
                                                   adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=lthf(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"lthf",buffer_2d(ni_start:ni_end, &
                                                   nj_start:nj_end),&
                                                   adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=fresh(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"fresh",buffer_2d(ni_start:ni_end, &
                                                     nj_start:nj_end),&
                                                     adios_err)
        error_sum = error_sum+adios_err

#ifdef COUP
!LPF20140119
          call exch_boundary(t_cpl,1)
          call exch_boundary(s_cpl,1)
          call exch_boundary(u_cpl,1)
          call exch_boundary(v_cpl,1)
          call exch_boundary(dhdx,1)
          call exch_boundary(dhdy,1)
          call exch_boundary(q,1)
!LPF20140119

        buffer_2d(:,:)=t_cpl(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"t_cpl",buffer_2d(ni_start:ni_end, &
                                                     nj_start:nj_end),&
                                                     adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=s_cpl(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"s_cpl",buffer_2d(ni_start:ni_end, &
                                                     nj_start:nj_end),&
                                                     adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=u_cpl(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"u_cpl",buffer_2d(ni_start:ni_end, &
                                                     nj_start:nj_end),&
                                                     adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=v_cpl(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"v_cpl",buffer_2d(ni_start:ni_end, &
                                                     nj_start:nj_end),&
                                                     adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=dhdx(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"dhdx",buffer_2d(ni_start:ni_end, &
                                                     nj_start:nj_end),&
                                                     adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=dhdy(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"dhdy",buffer_2d(ni_start:ni_end, &
                                                     nj_start:nj_end),&
                                                     adios_err)
        error_sum = error_sum+adios_err

        buffer_2d(:,:)=q(:,:)
        if(iy==0) then
        buffer_2d(:,1)=spval
        endif
        if(ix==0) then
        buffer_2d(1,:)=spval
        endif
        if(ix== (nx_proc-1)) then
        buffer_2d(imt,:)=spval
        endif
        call adios_write (adios_handle,"q",buffer_2d(ni_start:ni_end, &
                                                     nj_start:nj_end),&
                                                     adios_err)
        error_sum = error_sum+adios_err
#endif       
        if(error_sum .ne. 0) write(6,*) 'proc:',mytid,'error in writing content'

        call adios_close (adios_handle, adios_err)
        if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_close'

        if(mytid==0) write(*,*)'ok  fort.22'

        if ( mytid == 0 ) then
            open (17, file="rpointer.ocn.adios", form='formatted')
            write(17,'(a21)') fname
            close(17)
        end if
    end if

    return

end subroutine ssaveins_adios_noneoverlap

subroutine datachange(change,kk,mm,toru)
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use output_mod, only: spval

IMPLICIT NONE
integer :: kk,mm,toru
real(r4) :: change(imt,jmt,kk,mm)
       if (toru.eq.1) then  !t grid
           do m=1,mm
           do k=1,kk
            do i=1,imt
             do j=1,jmt
              if(vit(i,j,k)<0.5) then
               change(i,j,k,m)=spval
              endif
             enddo
            enddo
           end do
           end do
        else  !u grid
           do m=1,mm
           do k=1,kk
            do i=1,imt
             do j=1,jmt
              if(viv(i,j,k)<0.5) then
               change(i,j,k,m)=spval
              endif
             enddo
            enddo
           end do
           end do
        endif

   if(iy==0) then
    change(:,1,:,:)=spval 
   endif
   if(ix==0) then
    change(1,:,:,:)=spval
   endif
   if(ix==(nx_proc-1)) then
    change(imt,:,:,:)=spval 
   endif

 return

end subroutine datachange

subroutine ssavemon_adios_noneoverlap
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use output_mod
use dyn_mod
use tracer_mod
use cdf_mod
#ifdef LOWRES
use diag_mod
#endif
use msg_mod

use forc_mod
use work_mod
 IMPLICIT none
     
      logical :: hist_output,rest_output
      CHARACTER ( LEN =   4 ) :: ftail
     ! CHARACTER ( LEN =  24 ) :: fname
      CHARACTER ( LEN =  18 ) :: fname1
     ! CHARACTER ( LEN =  15 ) :: fname1
      CHARACTER ( LEN =   8 ) :: dd
      CHARACTER ( LEN =   10 ) :: tt
      CHARACTER ( LEN =   5 ) :: zz
      INTEGER(r4)             :: vv(8)
      INTEGER :: nwmf
!      real*4 vivtmp(imt,jmt,km) !LPF20131027
!      real*4 vittmp(imt,jmt,km) !LPF20131027

integer:: ni_offset,nj_offset,nk_offset, &
          ni_local,nj_local,nk1_local, &
          nk2_local,nk3_local,&
          ni_global,nj_global,nk1_global, &
          nk2_global,nk3_global, &
          ni_start,nj_start,nk_start, &
          ni_end,nj_end,nk1_end, &
          nk2_end, nk3_end, &
          adios_err,error_sum
integer*8:: adios_groupsize,adios_totalsize,adios_handle,group_comm
character (len=20) :: ttdd
integer :: i1,j1,k1
real*4 :: buffer_2d_tmp(imt,jmt)
real*4 :: buffer_tmp(imt,jmt,klv)
!---------------------------------------------------------------------
!     output monthly results
!---------------------------------------------------------------------
!    file name
      if (mon0-1 == 0 ) then
      nwmf = iyfm-1
      write (ftail,'(i4.4)') nwmf
      fname1(11:12)="12"
      else
      nwmf = iyfm
      write (ftail,'(i4.4)') nwmf
      write(fname1(11:12),'(i2.2)') mon0-1
      endif
 
!      write(fname1(11:12),'(i2.2)') iday
      fname1(1:5)='MMEAN'
      fname1(6:9)=ftail
      fname1(10:10)='-'
      !fname1(13:15)='.bp'

      fname1(13:13)='-'
      write(fname1(14:15),'(i2.2)') number_day
      fname1(16:18)='.bp'

!      if (iday==imd) then
      if (mod (iyfm,io_hist)==0 ) then
         hist_output=.true.
      else
        hist_output=.false.
      endif
!
      if (mod ((month-1),io_rest)==0 ) then
         rest_output=.true.
      else
         rest_output=.false.
      endif

       !write(*,*) 'hist_output,io_hist',hist_output,io_hist
       if(mytid==0) write(*,*) 'begin in ssavemon_adios_noneoverlap: ',fname1

   if (hist_output) then

      call adios_open (adios_handle,"ssavemon",fname1,"w",mpi_comm_ocn,adios_err)
        if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_open'

        ni_global= imt_global
        nj_global= jmt_global-jst_global+1
        nk1_global= km
        nk2_global= km+1
        nk3_global= klv

        if( nx_proc == 1 ) then
            ni_local = ni_global
            ni_offset= 0
            ni_start = 1
        else if( ix == 0 ) then
            ni_local = imt - num_overlap/2
            ni_offset= 0
            ni_start = 1
        else if( ix == (nx_proc-1) ) then
            ni_local = imt - num_overlap/2
            ni_offset= i_start(mytid+1)
            ni_start = num_overlap/2+1
        else
            ni_local = imt - num_overlap
            ni_offset= i_start(mytid+1)
            ni_start = num_overlap/2+1
        endif

        if( ny_proc == 1 ) then
            nj_local = nj_global
            nj_offset= 0
            nj_start = 1
        else if( iy == 0 ) then
            nj_local = jmt - num_overlap/2
            nj_offset= 0
            nj_start = 1
        else if( iy == (ny_proc-1) ) then
            nj_local = jmt_global - j_start(mytid+1)
            nj_offset= j_start(mytid+1)-jst_global+1
            nj_start = num_overlap/2+1
        else
            nj_local = jmt - num_overlap
            nj_offset= j_start(mytid+1)-jst_global+1
            nj_start = num_overlap/2+1
        endif

        nk1_local = km
        nk2_local = km+1
        nk3_local = klv
        nk_offset= 0
        nk_start= 1

        ni_end = ni_start+ni_local-1
        nj_end = nj_start+nj_local-1
        nk1_end = nk_start+nk1_local-1
        nk2_end = nk_start+nk2_local-1
        nk3_end = nk_start+nk3_local-1
        
        adios_groupsize= 4*13 + 8 + &
                         4*(ni_local+nj_local+2*km+1+1) + &
                         4*(14*ni_local*nj_local)+ &
                         4*(11*ni_local*nj_local*klv)
#if ( defined TIDEMIX )
       adios_groupsize= adios_groupsize+4*ni_local*nj_local*klv
#endif
#ifdef ISOOUT
        adios_groupsize= adios_groupsize+4*3*ni_local*nj_local*klv
#endif
        call adios_group_size (adios_handle,adios_groupsize, &
                               adios_totalsize,adios_err)
        if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_group_size'

        error_sum=0
        call adios_write(adios_handle,"ni_global",ni_global,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nj_global",nj_global,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk1_global",nk1_global,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk2_global",nk2_global,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk3_global",nk3_global,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"ni_offset",ni_offset,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nj_offset",nj_offset,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk_offset",nk_offset,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"ni_local",ni_local,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nj_local",nj_local,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk1_local",nk1_local,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk2_local",nk2_local,adios_err)
        error_sum = error_sum+adios_err
        call adios_write(adios_handle,"nk3_local",nk3_local,adios_err)
        error_sum = error_sum+adios_err

        if(error_sum .ne. 0) write(6,*) 'proc:',mytid,'error in writing dimensions'

        error_sum=0
        t0_cdf=month-1
        call adios_write (adios_handle,"lon",lon((ni_offset+1):(ni_offset+ni_local)),adios_err)
        error_sum = error_sum+adios_err
        call adios_write (adios_handle,"lat",lat((nj_offset+3):(nj_offset+3+nj_local)),adios_err)
        error_sum = error_sum+adios_err
        call adios_write (adios_handle,"lev",lev(nk_start:nk1_end),adios_err)
        error_sum = error_sum+adios_err
        call adios_write (adios_handle,"lev1",lev1(nk_start:nk2_end),adios_err)
         error_sum = error_sum+adios_err 
        call adios_write (adios_handle,"time",t0_cdf,adios_err)
         error_sum = error_sum+adios_err

         call adios_write (adios_handle,"spval",spval,adios_err)
         error_sum = error_sum+adios_err

        call date_and_time (dd,tt,zz,vv)
        ttdd=tt//'  '//dd
         call adios_write (adios_handle,"ttdd",ttdd,adios_err)

 
        call datachange(z0mon,1,1,1)
        call adios_write (adios_handle,"z0",z0mon(ni_start:ni_end, &
                                               nj_start:nj_end),&
                                               adios_err)
        error_sum = error_sum+adios_err
         call datachange(himon,1,1,1)
        call adios_write (adios_handle,"hi",himon(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err
      
        call datachange(hdmon,1,1,1)
        call adios_write (adios_handle,"hd",hdmon(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err
 
        DO j = 1,jmt
            DO i = 1,imt
               IF (vit (i,j,1) > .5) THEN
!                  icmon (i,j,1)= icmon (i,j,1)/Num_output
                   buffer_2d_tmp (i,j)= icmon (i,j,1)/nmonth (mon0)
               ELSE
                   buffer_2d_tmp (i,j)= spval
               END IF
            END DO
         END DO
        call datachange(buffer_2d_tmp,1,1,1) 
        call adios_write (adios_handle,"ic1",buffer_2d_tmp(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err
      
          DO j = 1,jmt
            DO i = 1,imt
               IF (vit (i,j,1) > .5) THEN
!                  icmon (i,j,2)= icmon (i,j,2)/Num_output
                   !buffer_2d_tmp (i,j)= icmon (i,j,2)/nmonth (mon0)
                   buffer_2d_tmp (i,j)= icmon (i,j,2)
               ELSE
                  buffer_2d_tmp (i,j)= spval
               END IF
            END DO
         END DO
         call datachange(buffer_2d_tmp,1,1,1)
         call adios_write (adios_handle,"ic2",buffer_2d_tmp(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err
        
         buffer_2d_tmp(:,:)=netmon(:,:,1)/ODZP(1)/OD0CP !W/m2
         call datachange(buffer_2d_tmp,1,1,1)
         call adios_write (adios_handle,"net1",buffer_2d_tmp(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err 
       
         buffer_2d_tmp(:,:)=netmon(:,:,2)*DZP(1)*1000./37.4*D0 !kg/m^2/s
         call datachange(buffer_2d_tmp,1,1,1)
         call adios_write (adios_handle,"net2",buffer_2d_tmp(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err

         call datachange(mldmon,1,1,1)
         call adios_write (adios_handle,"mld",mldmon(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err

!3D output
!         allocate(buffer3_r4(imt_global,jmt_global,km))
!         allocate(buffer3_r4msf(imt_global,jmt_global,km+1))
!         allocate(viv_tmp(imt_global,jmt_global,km)) !LPF20131026
!         allocate(vit_tmp(imt_global,jmt_global,km)) !LPF20131026
!         allocate(temp_tmp(imt_global,jmt_global,km)) !LPF20131030
!         allocate(salt_tmp(imt_global,jmt_global,km)) !LPF20131030
!         allocate(vel_tmp(imt_global,jmt_global,km)) !LPF20131030
!         allocate(veliso_tmp(imt_global,jmt_global,km)) !LPF20131030
!         allocate(ddy_tmp(imt_global,jmt_global,km,2)) !LPF20131030
!       vittmp=vit
!       vivtmp=viv
!       call local_to_global_4d(vittmp,vit_tmp,km,1,1)
!       where(vit_tmp(:,:,:)>10) vit_tmp=0.0
!       call local_to_global_4d(vivtmp,viv_tmp,km,1,0)
!       where(viv_tmp(:,:,:)>10) viv_tmp=0.0
!       call local_to_global_4d(tsmon,temp_tmp,km,1,1)
!       call local_to_global_4d(ssmon,salt_tmp,km,1,1)
!       call local_to_global_4d(vsmon,vel_tmp,km,1,0)

         call datachange(akmmon,klv,1,0)
         call adios_write (adios_handle,"akm",akmmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk2_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

         call datachange(aktmon,klv,1,1)
         call adios_write (adios_handle,"akt",aktmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

         call datachange(aksmon,klv,1,1)
         call adios_write (adios_handle,"aks",aksmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

          call datachange(akmbmon,klv,1,1)
         call adios_write (adios_handle,"akmb",akmbmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

         call datachange(aktbmon,klv,1,1)
         call adios_write (adios_handle,"aktb",aktbmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

         call datachange(aksbmon,klv,1,1)
         call adios_write (adios_handle,"aksb",aksbmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

#if ( defined TIDEMIX )
        call datachange(aktidebmon,klv,1,1)
         call adios_write (adios_handle,"aktide",aktidemon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err
#endif

        call datachange(tsmon,klv,1,1)
         call adios_write (adios_handle,"ts",tsmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

        call datachange(ssmon,klv,1,1)
         call adios_write (adios_handle,"ss",ssmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

        call datachange(wsmon,klv,1,1)
         call adios_write (adios_handle,"ws",wsmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

         call datachange(usmon,klv,1,0)
         call adios_write (adios_handle,"us",usmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

          call datachange(vsmon,klv,1,0)
         call adios_write (adios_handle,"vs",vsmon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err
 !         call local_to_global_4d(vsmon,buffer3_r4,klv,1,0)
 !         if(mytid==0) then
 !          if(diag_msf)  call msf(buffer3_r4,viv_tmp,psi)
 !          write(*,*)'after diag_msf'
 !         endif

         call datachange(sumon,1,1,0)
         call adios_write (adios_handle,"su",sumon(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err

         call datachange(svmon,1,1,0)
         call adios_write (adios_handle,"sv",svmon(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err

          call datachange(lthfmon,1,1,1)
         call adios_write (adios_handle,"lthf",lthfmon(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err

         call datachange(sshfmon,1,1,1)
         call adios_write (adios_handle,"sshf",sshfmon(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err

          call datachange(lwvmon,1,1,1)
         call adios_write (adios_handle,"lwv",lwvmon(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err

          call datachange(swvmon,1,1,1)
         call adios_write (adios_handle,"swv",swvmon(ni_start:ni_end, &
                                             nj_start:nj_end), &
                                             adios_err)
         error_sum = error_sum+adios_err

#ifdef ISOOUT
          call datachange(vetisomon,klv,1,1)
         call adios_write (adios_handle,"ustar",vetisomon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err

         do k=1,nk3_local
           do j=1,jmt
            do i=1,imt
             buffer_tmp(i,j,k) = (-1)* vntisomon(i,j,k)
            enddo
           enddo
         enddo
         call datachange(buffer_tmp,klv,1,1)
         call adios_write (adios_handle,"vstar",buffer_tmp(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err
 
          call datachange(vbtisomon,klv,1,1)
         call adios_write (adios_handle,"wstar",vbtisomon(ni_start:ni_end, &
                                             nj_start:nj_end, &
                                             nk_start:nk3_end),&
                                             adios_err)
         error_sum = error_sum+adios_err
#endif

!     if(diag_mth) then
!         if(mytid==0) write(*,*)'before diag_mth'
!LPF 20131223
!#ifdef ISOOUT
!       call diag_heat_transport(1,viv_tmp,vit_tmp,temp_tmp,salt_tmp,vel_tmp,veliso_tmp*(-1),ddy_tmp)
!       call diag_heat_transport(2,viv_tmp,vit_tmp,temp_tmp,salt_tmp,vel_tmp,veliso_tmp*(-1),ddy_tmp)
!#else
!        veliso_tmp=0.0
!       call diag_heat_transport(1,viv_tmp,vit_tmp,temp_tmp,salt_tmp,vel_tmp,veliso_tmp*(-1),ddy_tmp)
!       call diag_heat_transport(2,viv_tmp,vit_tmp,temp_tmp,salt_tmp,vel_tmp,veliso_tmp*(-1),ddy_tmp)
!#endif
!LPF 20131223
 !      call diag_heat_gyre_overturning(1,viv_tmp,vit_tmp,temp_tmp,salt_tmp,vel_tmp)
 !      call diag_heat_gyre_overturning(2,viv_tmp,vit_tmp,temp_tmp,salt_tmp,vel_tmp)
 !        if(mytid==0) write(*,*)'after diag_mth'
 !     endif
!LPF20131030
!LPF 20120728
 !     if(mytid==0) write(*,*)'before diag_msf'
!LPF20131024
!LPF20131024
 !     if(diag_bsf) then
 !        call barosf(usmon,bsf)
 !        if(mytid==0) write(*,*)'ok bsf'
!#ifdef ISOOUT
!         call barosf(vetisomon,bsf_iso)
!         if(mytid==0) write(*,*)'ok bsf_iso'
!#endif
!      endif
        
!         buffer3_r4msf(:,:,:)= 0.0
!        if (diag_msf) then
!          do i=1,3
!          do k=1,km+1
!          do j=1,jmt_global
             !if (psi(j,k,1)<10000.) then
                !t2z_cdf(j,k,1)=psi(j,k,1) !/(nmonth (mon0))/NSS
!                buffer3_r4msf(i,j,k)=psi(j,k,i)
!#ifdef ISOOUT
!                buffer3_r4msf(10+i,j,k)=psi_iso(j,k,i)
!#endif
             !end if
 !         end do
 !         end do
 !         end do

  !       call adios_write (adios_handle,"psi",buffer3_r4msf((ni_offset+ni_start):(ni_offset+ni_end), &
 !                                            (nj_offset+nj_start):(nj_offset+nj_end), &
 !                                            nk_start:nk2_end),&
 !                                            adios_err)
 !        error_sum = error_sum+adios_err
 !       endif !diag_msf
        
!       if (diag_bsf) then
!          do j=1,jmt_global
!          do i=1,imt_global
!                t2_cdf(i,j,1)=bsf(i,j) !/(nmonth (mon0))/NSS
!          end do
!          end do
!          call adios_write (adios_handle,"bsf",t2_cdf((ni_offset+ni_start):(ni_offset+ni_end), &
!                                             (nj_offset+nj_start):(nj_offset+nj_end),1),&
!                                             adios_err)
!         error_sum = error_sum+adios_err
!#ifdef ISOOUT
!          do j=1,jmt_global
!          do i=1,imt_global
!                t2_cdf(i,j,1)=bsf_iso(i,j) !/(nmonth (mon0))/NSS
!          end do
!          end do
!         call adios_write (adios_handle,"bsfiso",t2_cdf((ni_offset+ni_start):(ni_offset+ni_end), &
!                                             (nj_offset+nj_start):(nj_offset+nj_end),1),&
!                                             adios_err)
!         error_sum = error_sum+adios_err
!#endif
 !     endif !diag_bsf

 !      if (diag_mth) then
!!
!mth-temp
 !         do i=1,3
 !         do j=1,jmt_global
 !               t1_cdf(i,j,1)=mth(j,i,1)
 !         end do
 !         end do
!mth-salt
 !         do i=4,6
 !         do j=1,jmt_global
 !               t1_cdf(i,j,1)=mth(j,i-4+1,2)
 !         end do
 !         end do
!mth-advtemp
 !         do i=11,13
 !         do j=1,jmt_global
 !               t1_cdf(i,j,1)=mth_adv(j,i-11+1,1)
 !         end do
 !         end do
!mth-advsalt
 !         do i=14,16
 !         do j=1,jmt_global
 !               t1_cdf(i,j,1)=mth_adv(j,i-14+1,2)
 !         end do
 !         end do
!mth-isoadv temp
 !         do i=21,23
 !         do j=1,jmt_global
 !               t1_cdf(i,j,1)=mth_adv_iso(j,i-21+1,1)
 !         end do
 !         end do
!mth-isoadv salt
 !         do i=24,26
 !         do j=1,jmt_global
 !               t1_cdf(i,j,1)=mth_adv_iso(j,i-24+1,2)
 !         end do
 !         end do
!mth-diff temp
 !         do i=31,33
 !         do j=1,jmt_global
 !               t1_cdf(i,j,1)=mth_dif(j,i-31+1,1)
 !         end do
 !         end do
!mth-diff salt
  !        do i=34,36
  !        do j=1,jmt_global
  !              t1_cdf(i,j,1)=mth_dif(j,i-34+1,2)
  !        end do
  !        end do
!mth-gyre temp
  !        do i=41,43
  !        do j=1,jmt_global
  !              t1_cdf(i,j,1)=mth_gyre(j,i-41+1,1)
  !        end do
  !        end do
!mth-gyre salt
   !       do i=44,46
   !       do j=1,jmt_global
   !             t1_cdf(i,j,1)=mth_gyre(j,i-44+1,2)
   !       end do
   !       end do
!mth-ovrt temp
   !       do i=51,53
   !       do j=1,jmt_global
   !             t1_cdf(i,j,1)=mth_ovrt(j,i-51+1,1)
   !       end do
   !       end do
!mth-ovrt salt
   !       do i=54,56
   !       do j=1,jmt_global
  !!!! !             t1_cdf(i,j,1)=mth_ovrt(j,i-54+1,2)
   !       end do
   !       end do
   !      call adios_write (adios_handle,"mth",t1_cdf((ni_offset+ni_start):(ni_offset+ni_end), &
   !                                          (nj_offset+nj_start):(nj_offset+nj_end),1),&
   !                                          adios_err)
   !      error_sum = error_sum+adios_err
   !   endif !diag_mth

!#if (defined SMAG_OUT)
!$OMP PARALLEL DO PRIVATE (K,J,I)
!         DO k = 1,klv
!            DO j = 1,jmt_global
!               DO i = 1,imt
!                  IF (vit_global (i,j,k) > 0.5) THEN
!                     t3_cdf (i,j,k,1)= am3mon_io (i,j,k)/ (nmonth (mon0))/NSS
!                  ELSE
!                     t3_cdf (i,j,k,1)= spval
!                  END IF
!               END DO
!            END DO
!         END DO
!       call adios_write (adios_handle,"am",t3_cdf(ni_offset+ni_start:ni_offset+ni_end, &
!                                             nj_offset+nj_start:nj_offset+nj_end,&
!                                             nk_start:nk_end),&
!                                             adios_err)
!         error_sum = error_sum+adios_err
!#endif
        if(error_sum .ne. 0) write(6,*) 'proc:',mytid,'error in writing content'

        call adios_close (adios_handle, adios_err)
        if(adios_err .ne. 0) write(6,*) 'proc:',mytid,'error in adios_close'

      if(mytid==0) write(*,*)'ok  MMEAN.bp'

    !  deallocate (buffer3_r4,buffer3_r4msf)
    !  deallocate(viv_tmp) !LPF20131026
    !  deallocate(vit_tmp,vel_tmp,veliso_tmp) !LPF20131026
    !  deallocate(temp_tmp,salt_tmp) !LPF20131026
    !  deallocate(ddy_tmp) !LPF20131026

    endif !hist_output

     CALL mm00 (klv)
!lhl20120728      IF (mod ( (month -1),12) == 0)THEN
     IF (mod ( (month -1),io_rest) == 0)THEN
         CALL yy00
     END IF

        RETURN
end subroutine ssavemon_adios_noneoverlap

