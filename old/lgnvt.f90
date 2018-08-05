!
!  File name: lg.f90
!  Date:      2011/05/15 23:49
!  Author:    Lukas Vlcek
! 


program lg
    use rndmod
    implicit none
    logical :: linp
    character*80 :: str, fildmp, fnam
    integer*4 :: nx, ny, nz, ix, iy, iz, ixm, ixp, iym, iyp, izm, izp, itprn, ierr
    integer*4 :: jx, jy, jz, kx, ky, kz, js, jq, jn, ks, kq, kn, q, tq, i, iq, jj, kk
    integer*4 :: is, nc, it, im, xq, yq, ic, wx, wy, wz, iip, iin
    real*8 :: c, kTc, rw, dummy, wff, wss, wfs, frac1, h, b, Tr, rnn, dti, su, time, timax, tiprn
    real*8, dimension(:), allocatable :: dtq
    integer*4, dimension(:), allocatable :: nm
    integer*4, dimension(:,:,:), allocatable :: ni
    integer*4, dimension(:,:), allocatable :: grd
    integer*4, dimension(:,:,:,:), allocatable :: fi

    ! initialize random number generator 
    seed = 113314724325.0

    ! read simulation parameters
    read(*,*) str, wff, wss, wfs  ! relative fluid-wall/fluid-fluid interaction
    read(*,*) str, Tr  ! reduced temperature
    read(*,*) str, nx, ny, nz
    read(*,*) str, linp
    if (linp) then
        read(*,*) str, fildmp
    else
        read(*,*) str, frac1  ! fraction of phase 1
    end if
    read(*,*) str, timax, tiprn

    nc = 6
    c = real(nc) ! basic cubic lattice
    b = 0.22163/Tr

    allocate(fi(3,nx,ny,nz))
    fi = 0
    allocate(nm(0:nc), ni(3, nx*ny*nz, 0:nc), dtq(0:nc))
    nm = 0 ; ni = 0 ; dtq = 0.0

    ! intial configuration
    print *, 'initial configuration...' 
    if (linp) then
        open(1, file=fildmp, form='unformatted')
        read(1) fi
        close(1, status='keep')
    else
        iip = 0
        iin = 0
        do iz = 1, nz
            do iy = 1, ny
                do ix = 1, nx
                    if (rnd(dummy) < frac1) then
                        iip = iip + 1
                        fi(1, ix, iy, iz) = 1
                        if (iip > 0.5*nx*ny*nz) then
                            iin = iin + 1
                            fi(1, ix, iy, iz) = 0
                        end if
                    else
                        iin = iin + 1
                        fi(1, ix, iy, iz) = 0
                        if (iin > 0.5*nx*ny*nz) then
                            iip = iip + 1
                            fi(1, ix, iy, iz) = 1
                        end if
                    end if
                end do
            end do
        end do
    end if

    allocate(grd(3,6))
    grd = -1
    grd(1,1) = -2
    grd(1,2) =  0
    grd(2,3) = -2
    grd(2,4) =  0
    grd(3,5) = -2
    grd(3,6) =  0

    ! create neighbor-count lists
    print *, 'filling lattices and lists' 
    do iz = 1, nz
        izm = modulo(iz-2, nz) + 1
        izp = modulo(iz, nz) + 1
        do iy = 1, ny
            iym = modulo(iy-2, ny) + 1
            iyp = modulo(iy, ny) + 1
            do ix = 1, nx
                ixm = modulo(ix-2, nx) + 1
                ixp = modulo(ix, nx) + 1
                iq = fi(1, ix, iy, iz)
                is = 0
                if (fi(1, ixm, iy, iz) == iq) is = is + 1
                if (fi(1, ixp, iy, iz) == iq) is = is + 1
                if (fi(1, ix, iym, iz) == iq) is = is + 1
                if (fi(1, ix, iyp, iz) == iq) is = is + 1
                if (fi(1, ix, iy, izm) == iq) is = is + 1
                if (fi(1, ix, iy, izp) == iq) is = is + 1
                !is = fi(1,ix, iy, iz)*is + (1-fi(1,ix,iy,iz))*(nc - is)
                nm(is) = nm(is) + 1
                fi(2, ix, iy, iz) = is
                fi(3, ix, iy, iz) = nm(is)
                ni(1, nm(is), is) = ix
                ni(2, nm(is), is) = iy
                ni(3, nm(is), is) = iz
            end do
        end do
    end do


    do i = 0, nc-1
        dtq(i) = (1.0 - real(i)/c)*exp(-4.0*real(i)*b*wfs)
    end do
    dtq(nc) = 0.0

    print *, 'start ...'
    time = 0.0
    it = 0
    ic = 0
    do 
        ic = ic + 1
        dti = sum(dtq(0:5)*nm(0:5))

        rnn = rnd(dummy)*dti
        iq = 0
        su = 0.0
        do 
            su = su + dtq(iq)*nm(iq)
            if (su > rnn) exit
            iq = iq + 1
        end do
        if (iq > nc-1) iq = nc-1

        im = int(real(nm(iq))*rnd(dummy)) + 1
        if (im < 1) im = 1
        if (im > nm(iq)) im = nm(iq)
    
        ix = ni(1,im,iq)
        iy = ni(2,im,iq)
        iz = ni(3,im,iq)
        if (iq /= fi(2, ix, iy, iz) .or. im /= fi(3, ix, iy, iz)) then
            print *, ic
            print *,'iq', iq, fi(2, ix, iy, iz)
            print *,'im', im, fi(3, ix, iy, iz)
        end if
        is = fi(1, ix, iy, iz)
        !print *, 'is, iq, im', is, iq, im

        tq = int(rnd(dummy)*real(nc - iq)) + 1
        if (tq < 1) tq = 1
        if (tq > nc-iq) tq = nc - iq

        q = 0
        do jj = 1, nc
            jx = modulo(ix + grd(1, jj), nx) + 1
            jy = modulo(iy + grd(2, jj), ny) + 1
            jz = modulo(iz + grd(3, jj), nz) + 1
            js = fi(1, jx, jy, jz)
            jq = fi(2, jx, jy, jz)
            jn = fi(3, jx, jy, jz)

            xq = jq + 2*(is + js) - 4*is*js - 1

            if (is /= js) then          ! same sigma
                q = q + 1
                if (q == tq) then
                    fi(1, jx, jy, jz) = 1 - js
                    xq = nc - jq - 1

                    do kk = 1, nc
                        kx = modulo(jx + grd(1, kk), nx) + 1
                        ky = modulo(jy + grd(2, kk), ny) + 1
                        kz = modulo(jz + grd(3, kk), nz) + 1
                        ks = fi(1, kx, ky, kz) 
                        kq = fi(2, kx, ky, kz) 
                        kn = fi(3, kx, ky, kz) 

                        yq = kq + 2*(js + ks) - 4*js*ks - 1
                        if (ix == kx .and. iy == ky .and. iz == kz) then
                            fi(1, ix, iy, iz) = 1 - is
                            yq = nc - iq - 1
                        end if

                        if (yq /= kq) then
                            fi(2, kx, ky, kz) = yq
                            nm(yq) = nm(yq) + 1
                            fi(3, kx, ky, kz) = nm(yq) 
                            ni(1, nm(yq), yq) = kx ! end of stack
                            ni(2, nm(yq), yq) = ky
                            ni(3, nm(yq), yq) = kz
                            if (kn /= nm(kq)) then ! fill the hole
                                wx = ni(1, nm(kq), kq)
                                wy = ni(2, nm(kq), kq)
                                wz = ni(3, nm(kq), kq)
                                ni(1, kn, kq) = wx
                                ni(2, kn, kq) = wy
                                ni(3, kn, kq) = wz
                                fi(3, wx, wy, wz) = kn
                            end if
                            nm(kq) = nm(kq) - 1
                        end if
                    end do
                end if
            end if

            jq = fi(2, jx, jy, jz)
            jn = fi(3, jx, jy, jz)
            if (xq /= jq) then
                fi(2, jx, jy, jz) = xq
                nm(xq) = nm(xq) + 1
                fi(3, jx, jy, jz) = nm(xq) 
                ni(1, nm(xq), xq) = jx ! end of stack
                ni(2, nm(xq), xq) = jy
                ni(3, nm(xq), xq) = jz
                if (jn /= nm(jq)) then ! fill the hole
                    wx = ni(1, nm(jq), jq)
                    wy = ni(2, nm(jq), jq)
                    wz = ni(3, nm(jq), jq)
                    ni(1, jn, jq) = wx
                    ni(2, jn, jq) = wy
                    ni(3, jn, jq) = wz
                    fi(3, wx, wy, wz) = jn
                end if
                nm(jq) = nm(jq) - 1
            end if
        end do

        time = time + 1/dti

        if (time-tiprn*real(it) > tiprn) then
            print *,'time',time,nm(0:6),sum(fi(1,:,:,:))/real(nx*ny*nz)
            it = int(time/tiprn)
            if (mod(it, 100) == 0) then
                write(fnam, '(i3)') it/100
                fnam = adjustl(fnam)
                ierr = len_trim(fnam)
                fnam(ierr+1:ierr+4) = '.dmp'
                open(2, file=fnam, form='unformatted', status='unknown')
                write(2) fi(1,:,:,:)
                close(2, status='keep')
            end if
        end if
        if (time > timax) exit
    end do

    stop
end program lg
