module mod_main

implicit none

contains


subroutine sus_iter(Gs,Ps,knotsm,knotsn,mmat0,nmat0,xmat0,pmat0,vmat0,wvecd,bet,kap,invlam,sig,rstar,checkMult, &
    mmat1,nmat1,xmat1,pmat1,rmat1,vmat1,zmat1,idmat)

    real(8), intent(in)  :: Gs(:), Ps(:,:), knotsm(:), knotsn(:)
    real(8), intent(in)  :: mmat0(:,:,:), nmat0(:,:,:), xmat0(:,:,:), pmat0(:,:,:), vmat0(:,:,:), wvecd(:)
    real(8), intent(in)  :: bet, kap, invlam, sig, rstar
    integer, intent(in)  :: checkMult
    real(8), intent(out) :: mmat1(:,:,:), nmat1(:,:,:), xmat1(:,:,:), pmat1(:,:,:), rmat1(:,:,:), &
    vmat1(:,:,:), zmat1(:,:,:), idmat(:,:,:)
    integer ns, nm, nn
    integer is, js, im, in, id, flagSUS, flagZLB
    real(8) snow, wstar, mpast, npast, m1, n1, x1, p1, r1
    real(8) ex, edxm, edxn, ep, edpm, edpn, xx1(1), xx2(2)
    ! real(8) vmat0(ns,nm,nn), vmat1(ns,nm,nn)
    ! real(8) zmat0(ns,nm,nn), zmat1(ns,nm,nn), idmat(ns,nm,nn), wvecd(ns)
    real(8) M, N, ev, edvm, edvn, v1, z1, xx3(1), xx4(2)
    real(8) ml, pl, xl, rl, mr, pr, xr, rr
    real(8), allocatable :: xcond(:,:), pcond(:,:), vcond(:,:)

    ns = size(Gs,1)
    nm = size(knotsm,1)
    nn = size(knotsn,1)
    allocate(xcond(nm,nn),pcond(nm,nn),vcond(nm,nn))


!    checkMult = 0

    do is = 1,ns

        snow = Gs(is)
        wstar = wvecd(is)

        xcond = 0.0d0
        pcond = 0.0d0
        vcond = 0.0d0

        do js = 1,ns

            xcond = xcond + Ps(is,js)*reshape(xmat0(js,:,:),(/nm,nn/))
            pcond = pcond + Ps(is,js)*reshape(pmat0(js,:,:),(/nm,nn/))
            vcond = vcond + Ps(is,js)*reshape(vmat0(js,:,:),(/nm,nn/))

        end do

        do im = 1,nm

            mpast = knotsm(im)

            do in = 1,nn

                npast = knotsn(in)
                M = -(1.0d0/bet/sig)*npast + mpast
                N = (1.0d0/bet)*npast

                ! initialize the flags
                flagSUS = 0
                flagZLB = 0

                ! check the sustainability constraint
!                if (vmat0(is,im,in)<=wstar) then
                if (vmat0(is,im,in)<=wstar .and. is==1) then ! 11/21/16: is==1 is not needed after fixing the FOCs?
                    flagSUS = 1
                end if

                if (flagSUS==0) then

                    id = 1
                    xx1(1) = mmat0(is,im,in)
                    xx1 = nra2(xx1,1d-6,50,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,0.0d0, &
                    bet,kap,invlam,sig,rstar)
                    m1 = xx1(1)
                    n1 = 0.0d0

                    call intf2(knotsm,knotsn,xcond,m1,n1,ex,edxm,edxn)
                    call intf2(knotsm,knotsn,pcond,m1,n1,ep,edpm,edpn)
                    x1 = (1.0d0/bet)*npast - kap*m1
                    p1 = kap*x1 + bet*ep ! if kap=0 and ep=0, p0=0
                    r1 = sig*(ex-x1) + ep + snow ! x1 = ex - (1/sig)*(r1-ep-snow)

                    call intf2(knotsm,knotsn,vcond,m1,n1,ev,edvm,edvn)
                    v1 = -x1**2 - invlam*p1**2 + bet*ev
                    z1 = 1.0d0

                    ! check the ZLB
                    if (r1<=-rstar) then

                        flagZLB = 1

                    end if

                    if (flagZLB==1) then

                        id = 2
                        xx2(1) = mmat0(is,im,in)
                        xx2(2) = nmat0(is,im,in)
                        xx2 = nra2(xx2,1d-6,50,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,0.0d0,bet,kap,invlam,sig,rstar)
                        m1 = xx2(1)
                        n1 = xx2(2)

                        call intf2(knotsm,knotsn,pcond,m1,n1,ep,edpm,edpn)
                        x1 = -n1 + (1.0d0/bet)*npast - kap*m1
                        p1 = kap*x1 + bet*ep ! if kap=0 and ep=0, p0=0
                        r1 = -rstar

                        call intf2(knotsm,knotsn,vcond,m1,n1,ev,edvm,edvn)
                        v1 = -x1**2 - invlam*p1**2 + bet*ev
                        z1 = 1.0d0

                    end if

                else

                    id = 3
                    xx3(1) = mmat0(is,im,in)
                    xx3 = nra2(xx3,1d-6,50,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar)
                    m1 = xx3(1)
                    n1 = 0.0d0

                    call intf2(knotsm,knotsn,xcond,m1,0.0d0,ex,edxm,edxn)
                    call intf2(knotsm,knotsn,pcond,m1,0.0d0,ep,edpm,edpn)
!                    p1 = (-kap*M*n1 + kap*(N-kap*M)*m1 + bet*M*ep)/(M+kap*invlam*N)
                    p1 = (kap*(N-kap*M)*m1 + bet*M*ep)/(M+kap*invlam*N)
                    x1 = (N/M-kap)*m1 - invlam*N/M*p1
                    r1 = sig*(ex-x1) + ep + snow ! x1 = ex - (1/sig)*(r1-ep-snow)

                    v1 = wstar
!                    z1 = bet*(x1+n1+kap*m1)/npast
                    z1 = bet*(x1+kap*m1)/max(npast,0.0001d0)

                    ! check the ZLB: a possibility of multiple solution ...
                    if (r1<=-rstar) then

                        flagZLB = 1;

                    elseif (checkMult==1) then ! check the other solutions

                        ml = m1
                        ml = gss(knotsm(1),m1+1d-3,1d-6,50,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar)
!                        ml = gss(knotsm(1)-10.0d0,m1+1d-3,1d-6,50,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar)
                        call intf2(knotsm,knotsn,xcond,ml,0.0d0,ex,edxm,edxn)
                        call intf2(knotsm,knotsn,pcond,ml,0.0d0,ep,edpm,edpn)
                        pl = (-kap*M*n1 + kap*(N-kap*M)*ml + bet*M*ep)/(M+kap*invlam*N)
                        xl = (N/M-kap)*ml - invlam*N/M*pl
                        rl = sig*(ex-xl) + ep + snow ! x1 = ex - (1/sig)*(r1-ep-snow)

                        mr = m1
                        mr = gss(m1-1d-3,knotsm(nm),1d-6,50,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar)
                        call intf2(knotsm,knotsn,xcond,mr,0.0d0,ex,edxm,edxn)
                        call intf2(knotsm,knotsn,pcond,mr,0.0d0,ep,edpm,edpn)
                        pr = (-kap*M*n1 + kap*(N-kap*M)*mr + bet*M*ep)/(M+kap*invlam*N)
                        xr = (N/M-kap)*mr - invlam*N/M*pr
                        rr = sig*(ex-xr) + ep + snow ! x1 = ex - (1/sig)*(r1-ep-snow)

                        if (rl<=-rstar .or. rr<=-rstar) then

                            flagZLB = 1

                        end if

                    end if

                    if (flagZLB==1) then

                        id = 4
                        xx4(1) = mmat0(is,im,in)
                        xx4(2) = nmat0(is,im,in)
                        xx4 = nra2(xx4,1d-6,50,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar)
                        m1 = xx4(1)
                        n1 = xx4(2)

                        call intf2(knotsm,knotsn,pcond,m1,n1,ep,edpm,edpn)
                        p1 = (-kap*M*n1 + kap*(N-kap*M)*m1 + bet*M*ep)/(M+kap*invlam*N)
                        x1 = -n1 + (N/M-kap)*m1 - invlam*N/M*p1
                        r1 = -rstar

                        v1 = wstar
                        z1 = bet*(x1+n1+kap*m1)/max(npast,0.0001d0)

                    end if

                end if

                mmat1(is,im,in) = m1
                nmat1(is,im,in) = n1
                xmat1(is,im,in) = x1
                pmat1(is,im,in) = p1
                rmat1(is,im,in) = r1
                vmat1(is,im,in) = v1
                zmat1(is,im,in) = z1
                idmat(is,im,in) = id

            end do

        end do

    end do


end subroutine sus_iter


subroutine comm_iter(Gs,Ps,knotsm,knotsn,mmat0,nmat0,xmat0,pmat0,vmat0,bet,kap,invlam,sig,rstar,mmat1,nmat1,xmat1,pmat1,rmat1,vmat1)

    implicit none

    real(8), intent(in)  :: Gs(:), Ps(:,:), knotsm(:), knotsn(:)
    real(8), intent(in)  :: mmat0(:,:,:), nmat0(:,:,:), xmat0(:,:,:), pmat0(:,:,:), vmat0(:,:,:)
    real(8), intent(in)  :: bet, kap, invlam, sig, rstar
    real(8), intent(out) :: mmat1(:,:,:), nmat1(:,:,:), xmat1(:,:,:), pmat1(:,:,:), rmat1(:,:,:), vmat1(:,:,:)
    integer ns, nm, nn
    ! real(8) nmin, nmax
    integer is, js, im, in, id
    real(8) snow, mpast, npast, m1, n1, x1, p1, r1, v1
    real(8) ex, edxm, edxn, ep, edpm, edpn, ev, edvm, edvn, xx1(1), xx2(2)
    real(8), allocatable :: xcond(:,:), pcond(:,:), vcond(:,:)

    ns = size(Gs,1)
    nm = size(knotsm,1)
    nn = size(knotsn,1)
    allocate(xcond(nm,nn),pcond(nm,nn),vcond(nm,nn))


    do is = 1,ns

        snow = Gs(is)

        xcond = 0.0d0
        pcond = 0.0d0
        vcond = 0.0d0

        do js = 1,ns

            xcond = xcond + Ps(is,js)*reshape(xmat0(js,:,:),(/nm,nn/))
            pcond = pcond + Ps(is,js)*reshape(pmat0(js,:,:),(/nm,nn/))
            vcond = vcond + Ps(is,js)*reshape(vmat0(js,:,:),(/nm,nn/))

        end do

        do im = 1,nm

            mpast = knotsm(im)

            do in = 1,nn

                npast = knotsn(in)

                id = 1
                xx1(1) = mmat0(is,im,in)
                xx1 = nra2(xx1,1d-6,50,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,0.0d0,bet,kap,invlam,sig,rstar)
                m1 = xx1(1)
                n1 = 0.0d0

                call intf2(knotsm,knotsn,xcond,m1,n1,ex,edxm,edxn)
                call intf2(knotsm,knotsn,pcond,m1,n1,ep,edpm,edpn)
                x1 = (1.0d0/bet)*npast - kap*m1
                p1 = kap*x1 + bet*ep ! if kap=0 and ep=0, p0=0
                r1 = sig*(ex-x1) + ep + snow ! x1 = ex - (1/sig)*(r1-ep-snow)

                if (r1<=-rstar) then

                    id = 2
                    xx2(1) = mmat0(is,im,in)
                    xx2(2) = nmat0(is,im,in)
                    xx2 = nra2(xx2,1d-6,50,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,0.0d0,bet,kap,invlam,sig,rstar);
                    m1 = xx2(1)
                    n1 = xx2(2)

                    call intf2(knotsm,knotsn,pcond,m1,n1,ep,edpm,edpn)
                    x1 = -n1 + (1.0d0/bet)*npast - kap*m1
                    p1 = kap*x1 + bet*ep ! if kap=0 and ep=0, p0=0
                    r1 = -rstar

                end if

                call intf2(knotsm,knotsn,vcond,m1,n1,ev,edvm,edvn)
                v1 = -x1**2 - invlam*p1**2 + bet*ev

                mmat1(is,im,in) = m1
                nmat1(is,im,in) = n1
                xmat1(is,im,in) = x1
                pmat1(is,im,in) = p1
                rmat1(is,im,in) = r1
                vmat1(is,im,in) = v1

            end do

        end do

    end do

end subroutine comm_iter


function gss(xmin,xmax,crit,maxiter,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar) result(x)

    real(8), intent(in)  :: xmin, xmax, crit
    integer, intent(in)  :: maxiter, id
    real(8), intent(in)  :: snow, mpast, npast, knotsm(:), knotsn(:), xcond(:,:), pcond(:,:), vcond(:,:), wstar
    real(8), intent(in)  :: bet, kap, invlam, sig, rstar

    real(8) rg, a, b, c, d, x, fc, fd, z
    real(8) diff
    integer iter


    rg = (3.0d0-sqrt(5.0d0))/2.0d0

    a = xmin
    b = xmax
    c = a + rg*(b-a)
    d = a + (1.0d0-rg)*(b-a)

    fc = focz2_gss(c,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar)
    fd = focz2_gss(d,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar)

    diff = 1d+4
    iter = 0

    x = c

    do while (diff>crit .and. iter<maxiter)

        if (fc>=fd) then

            z = c + (1.0d0-rg)*(b-c)
            a = c
            c = d
            fc = fd
            d = z
            fd = focz2_gss(d,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar)

        else

            z = a + rg*(d-a)
            b = d
            d = c
            fd = fc
            c = z
            fc = focz2_gss(c,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar)

        end if

        diff = d-c
        iter = iter+1

    end do

    x = c


end function gss


function focz2_gss(m0,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar) result(f)

    real(8), intent(in)  :: m0, snow, mpast, npast, knotsm(:), knotsn(:)
    real(8), intent(in)  :: xcond(:,:), pcond(:,:), vcond(:,:), wstar, bet, kap, invlam, sig, rstar
    real(8) f
    real(8) M, N, ex, edxm, edxn, ep, edpm, edpn, ev, edvm, edvn, x0, p0, dx0, dp0

    ! m0 = xx3(1)

    M = -(1.0d0/bet/sig)*npast + mpast
    N = (1.0d0/bet)*npast

    ! call intf2(knotsm,knotsn,xcond,m0,0.0d0,ex,edxm,edxn)
    call intf2(knotsm,knotsn,pcond,m0,0.0d0,ep,edpm,edpn)
    call intf2(knotsm,knotsn,vcond,m0,0.0d0,ev,edvm,edvn)

    p0 = (kap*(N-kap*M)*m0 + bet*M*ep)/(M+kap*invlam*N)
    dp0 = (kap*(N-kap*M) + bet*M*edpm)/(M+kap*invlam*N)
    x0 = (N/M-kap)*m0 - invlam*N/M*p0
    dx0 = (N/M-kap) - invlam*N/M*dp0

    f = -wstar -(x0**2+invlam*p0**2) + bet*ev
    ! df(1,1) = -2.0d0*x0*dx0 - invlam*2.0d0*p0*dp0 + bet*edvm

    ! for gss
    f = abs(f)

end function focz2_gss


function nra2(xinit,crit,maxiter,id,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar) result(x0)

    real(8), intent(in)  :: xinit(:), crit
    integer, intent(in)  :: maxiter, id
    real(8), intent(in)  :: snow, mpast, npast, knotsm(:), knotsn(:), xcond(:,:), pcond(:,:), vcond(:,:), wstar
    real(8), intent(in)  :: bet, kap, invlam, sig, rstar

    real(8) diff
    real(8), allocatable :: x0(:), dx(:), x1(:), f(:), df(:,:)
    integer iter, ix, nx


    nx = size(xinit,1)
    allocate(x0(nx),dx(nx),x1(nx),f(nx),df(nx,nx))
	x0 = xinit
	! xlow  = xmin
	! xhigh = xmax

	diff = 1d+4;
	iter = 0;

	! Newton-Rhapson
	do while (diff>crit .and. iter<maxiter)

		if (id==1) then
			call foc2(x0,snow,mpast,npast,knotsm,knotsn,xcond,pcond,bet,kap,invlam,sig,rstar,f,df)
        elseif (id==2) then
			call focr2(x0,snow,mpast,npast,knotsm,knotsn,xcond,pcond,bet,kap,invlam,sig,rstar,f,df);
        elseif (id==3) then
			call focz2(x0,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar,f,df);
        elseif (id==4) then
			call focrz2(x0,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar,f,df);
		end if

		if (mycond(df)<1d-12) then
            do ix = 1,nx
                df(ix,ix) = df(ix,ix) + 1d-6
            end do
        end if
	    dx = -matmul(myinv(df),f)
	    x1 = x0 + dx

		diff = maxval(abs(x1-x0),1)
	    iter = iter + 1
	    x0 = x1

	end do

	! // // check bounds
	! // if (x0<xmin) {
	! //     x0 = xmin;
	! // }
	! // else if (x0>xmax) {
	! //     x0 = xmax;
	! // }

	! return x0;


end function nra2


function mycond(A) result(cond)

    real(8), intent(in) :: A(:,:)
    real(8), allocatable :: invA(:,:)
    integer n
    real(8) normA, norminvA, cond

    n = size(A,1)
    allocate(invA(n,n))

    invA = myinv(A)

    if (n==1) then
        normA = sum(abs(A))
        norminvA = sum(abs(invA))
    else
        normA = maxval(sum(abs(A),1),1)
        norminvA = maxval(sum(abs(invA),1),1)
    end if

    cond = 1.0d0/normA/norminvA

end function mycond


function myinv(X) result(invX)

    real(8), intent(in) :: X(:,:)
    real(8), allocatable :: invX(:,:)
    real(8) a, b, c, d
    integer n

    n = size(X,1)
    allocate(invX(n,n))

    if (n==1) then
        invX(1,1) = 1.0d0/X(1,1)
    elseif (n==2) then
        a = X(1,1)
        b = X(1,2)
        c = X(2,1)
        d = X(2,2)
        invX(1,:) = (/d, -b/)
        invX(2,:) = (/-c, a/)
        invX = invX/(a*d-b*c)
    end if

end function myinv


subroutine foc2(xx1,snow,mpast,npast,knotsm,knotsn,xcond,pcond,bet,kap,invlam,sig,rstar,f,df)

    real(8), intent(in)  :: xx1(:), snow, mpast, npast, knotsm(:), knotsn(:), xcond(:,:), pcond(:,:), bet, kap, invlam, sig, rstar
    real(8), intent(out) :: f(:), df(:,:)
    real(8) m0, ex, edxm, edxn, ep, edpm, edpn, x0, p0, dx0, dp0

    m0 = xx1(1)

    call intf2(knotsm,knotsn,xcond,m0,0.0d0,ex,edxm,edxn)
    call intf2(knotsm,knotsn,pcond,m0,0.0d0,ep,edpm,edpn)

    x0 = (1.0d0/bet)*npast - kap*m0
    p0 = kap*x0 + bet*ep ! if kap=0 and ep=0, p0=0
    dx0 = -kap
    dp0 = kap*dx0 + bet*edpm

    f(1) = -invlam*p0 + (1.0d0/bet/sig)*npast + m0 - mpast ! if invlam=0, m0=mpast-(1/bet/sig)*npast
    df(1,1) = -invlam*dp0 + 1.0d0

end subroutine foc2


subroutine focr2(xx2,snow,mpast,npast,knotsm,knotsn,xcond,pcond,bet,kap,invlam,sig,rstar,f,df)

    real(8), intent(in)  :: xx2(:), snow, mpast, npast, knotsm(:), knotsn(:), xcond(:,:), pcond(:,:), bet, kap, invlam, sig, rstar
    real(8), intent(out) :: f(:), df(:,:)
    real(8) m0, n0, ex, edxm, edxn, ep, edpm, edpn, x0, dxm, dxn, p0, dpm, dpn, r0

    m0 = xx2(1)
    n0 = xx2(2)

    call intf2(knotsm,knotsn,xcond,m0,n0,ex,edxm,edxn)
    call intf2(knotsm,knotsn,pcond,m0,n0,ep,edpm,edpn)

    x0 = -n0 + (1.0d0/bet)*npast - kap*m0
    dxm = -kap
    dxn = -1.0d0
    p0 = kap*x0 + bet*ep ! if kap=0 and ep=0, p0=0
    dpm = kap*dxm + bet*edpm
    dpn = kap*dxn + bet*edpn
    r0 = -rstar

    f(1) = -invlam*p0 + (1.0d0/bet/sig)*npast + m0 - mpast ! if invlam=0, m0=mpast-(1/bet/sig)*npast
    f(2) = -x0 + ex - (1.0d0/sig)*(r0-ep-snow)
    df(1,1) = -invlam*dpm + 1.0d0
    df(1,2) = -invlam*dpn
    df(2,1) = -dxm + edxm - (1.0d0/sig)*(-edpm)
    df(2,2) = -dxn + edxn - (1.0d0/sig)*(-edpn)

end subroutine focr2


subroutine focz2(xx3,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar,f,df)

    real(8), intent(in)  :: xx3(:), snow, mpast, npast, knotsm(:), knotsn(:)
    real(8), intent(in)  :: xcond(:,:), pcond(:,:), vcond(:,:), wstar, bet, kap, invlam, sig, rstar
    real(8), intent(out) :: f(:), df(:,:)
    real(8) m0, M, N, ex, edxm, edxn, ep, edpm, edpn, ev, edvm, edvn, x0, p0, dx0, dp0

    m0 = xx3(1)

    M = -(1.0d0/bet/sig)*npast + mpast
    N = (1.0d0/bet)*npast

    ! call intf2(knotsm,knotsn,xcond,m0,0.0d0,ex,edxm,edxn)
    call intf2(knotsm,knotsn,pcond,m0,0.0d0,ep,edpm,edpn)
    call intf2(knotsm,knotsn,vcond,m0,0.0d0,ev,edvm,edvn)

    p0 = (kap*(N-kap*M)*m0 + bet*M*ep)/(M+kap*invlam*N)
    dp0 = (kap*(N-kap*M) + bet*M*edpm)/(M+kap*invlam*N)
    x0 = (N/M-kap)*m0 - invlam*N/M*p0
    dx0 = (N/M-kap) - invlam*N/M*dp0

    f(1) = -wstar -(x0**2+invlam*p0**2) + bet*ev
    df(1,1) = -2.0d0*x0*dx0 - invlam*2.0d0*p0*dp0 + bet*edvm

end subroutine focz2


subroutine focrz2(xx4,snow,mpast,npast,knotsm,knotsn,xcond,pcond,vcond,wstar,bet,kap,invlam,sig,rstar,f,df)

    real(8), intent(in)  :: xx4(:), snow, mpast, npast, knotsm(:), knotsn(:)
    real(8), intent(in)  :: xcond(:,:), pcond(:,:), vcond(:,:), wstar, bet, kap, invlam, sig, rstar
    real(8), intent(out) :: f(:), df(:,:)
    real(8) m0, n0, M, N, ex, edxm, edxn, ep, edpm, edpn, ev, edvm, edvn
    real(8) x0, dxm, dxn, p0, dpm, dpn, r0

    m0 = xx4(1)
    n0 = xx4(2)

    M = -(1.0d0/bet/sig)*npast + mpast
    N = (1.0d0/bet)*npast

    call intf2(knotsm,knotsn,xcond,m0,n0,ex,edxm,edxn)
    call intf2(knotsm,knotsn,pcond,m0,n0,ep,edpm,edpn)
    call intf2(knotsm,knotsn,vcond,m0,n0,ev,edvm,edvn)

    p0 = (-kap*M*n0 + kap*(N-kap*M)*m0 + bet*M*ep)/(M+kap*invlam*N)
    x0 = -n0 + (N/M-kap)*m0 - invlam*N/M*p0 ! fixed on 11/22/16
    dpm = (kap*(N-kap*M) + bet*M*edpm)/(M+kap*invlam*N)
    dpn = (-kap*M + bet*M*edpn)/(M+kap*invlam*N)
    dxm = (N/M-kap) - invlam*N/M*dpm
    dxn = -1 - invlam*N/M*dpn
    r0 = -rstar

    f(1) = -wstar -(x0**2+invlam*p0**2) + bet*ev
    f(2) = -x0 + ex - (1.0d0/sig)*(r0-ep-snow)
    df(1,1) = -2.0d0*x0*dxm - invlam*2.0d0*p0*dpm + bet*edvm
    df(1,2) = -2.0d0*x0*dxn - invlam*2.0d0*p0*dpn + bet*edvn
    df(2,1) = -dxm + edxm - (1.0d0/sig)*(-edpm)
    df(2,2) = -dxn + edxn - (1.0d0/sig)*(-edpn)

end subroutine focrz2


subroutine intf2(xgrid,ygrid,zmat,x,y,z,dzx,dzy)

    real(8), intent(in)  :: xgrid(:), ygrid(:), zmat(:,:), x, y
    real(8), intent(out) :: z, dzx, dzy
    integer nx, ix, jx, ny, iy, jy
    real(8) etax, etay, z1, z2, detax, detay, dz1y, dz2y


	! 2-dim linear interpolation
    nx = size(xgrid,1)
    ny = size(ygrid,1)

    ! for x
	ix = 1
	if (x<=xgrid(1)) then
	    ix = 1
	elseif (x>=xgrid(nx)) then
	    ix = nx-1
	else
	    jx = 2
	    do while (jx<=nx)

	        if (x<xgrid(jx)) then
	            ix = jx-1
	            jx = nx
	        end if

	    jx = jx+1
        end do
	end if

    ix = min(max(1,ix),nx-1)

    ! for y
	iy = 1
	if (y<=ygrid(1)) then
	    iy = 1
	elseif (y>=ygrid(ny)) then
	    iy = nx-1
	else
	    jy = 2
	    do while (jy<=ny)

	        if (y<ygrid(jy)) then
	            iy = jy-1
	            jy = ny
	        end if

	    jy = jy+1
        end do
	end if

    iy = min(max(1,iy),ny-1)

	etax = (x-xgrid(ix))/(xgrid(ix+1)-xgrid(ix))
    etay = (y-ygrid(iy))/(ygrid(iy+1)-ygrid(iy))
    z1 = (1.0d0-etay)*zmat(ix,iy)   + etay*zmat(ix,iy+1)
    z2 = (1.0d0-etay)*zmat(ix+1,iy) + etay*zmat(ix+1,iy+1)
    z  = (1.0d0-etax)*z1 + etax*z2

    detax = 1.0d0/(xgrid(ix+1)-xgrid(ix))
    dzx = -detax*z1 + detax*z2

    detay = 1.0d0/(ygrid(iy+1)-ygrid(iy))
    dz1y = -detay*zmat(ix,iy) + detay*zmat(ix,iy+1)
    dz2y = -detay*zmat(ix+1,iy) + detay*zmat(ix+1,iy+1)
    dzy = (1.0d0-etax)*dz1y + etax*dz2y

end subroutine intf2


function linspace(a,b,m) result(x)

    integer, intent(in) :: m
    real(8), intent(in) :: a, b
    real(8), dimension(m) :: x
    integer i

    x(1) = a

    do i = 2,m

        x(i) = x(i-1) + (b-a)/(m-1)

    end do

end function linspace


function logspace(a,b,m) result(x)

    integer, intent(in) :: m
    real(8), intent(in) :: a, b
    real(8), dimension(m) :: x
    integer i

    x = linspace(log(10.0d0**a),log(10.0d0**b),m)
    x = exp(x)

end function logspace


end module mod_main
