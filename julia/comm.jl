# Optimal commitment policy
using NLsolve
using Interpolations

function main()

    damp = 1.0

    bet = 0.9925
    kap = 0.01
    phi = 5.0
    rstar = (1/bet-1)*100
    zeta = 1.0
    invlam = 1/0.003
    sig = 1.0

    pH = 0.0
    pL = 0.8
    # NOTE: i0 is deviation from the steady state
    sH = 0.0 #rstar
    sL = -2.5 #rstar-2.5

    # policy function iteration
    ns = 2
    nm = 21
    nn = 21
    mmax = 0.0
    mmin = -200.0
    nmax = 50.0
    nmin = 0.0

    Gs = [sH sL]'
    Ps = [1-pH pH; 1-pL pL]

    knotsm = range(mmin,stop=mmax,length=nm)
    knotsn = exp10.(range(log(1.0)/log(10.0),stop=log(nmax-nmin+1.0)/log(10.0),length=nn))
    knotsn = knotsn .+ (nmin - 1.0)

    xvec0, pvec0, rvec0, vvec0 = calcmpe(Gs,Ps,bet,kap,invlam,sig,rstar)
    println([xvec0[2] pvec0[2]])

    ymat0 = zeros(Float64,nm,nn,ns)
    pmat0 = zeros(Float64,nm,nn,ns)
    imat0 = zeros(Float64,nm,nn,ns)
    mmat0 = zeros(Float64,nm,nn,ns)
    nmat0 = zeros(Float64,nm,nn,ns)
    vmat0 = zeros(Float64,nm,nn,ns)
    ymat1 = copy(ymat0)
    pmat1 = copy(pmat0)
    imat1 = copy(imat0)
    mmat1 = copy(mmat0)
    nmat1 = copy(nmat0)
    vmat1 = copy(vmat0)

    crit = 1e-5
    diff = 1e+4
    iter = 0

    while (diff>crit)

        @inbounds for is in 1:ns

            s0 = Gs[is]

            eyvec = zeros(Float64,nm,nn)
            epvec = zeros(Float64,nm,nn)
            evvec = zeros(Float64,nm,nn)

            @inbounds for js in 1:ns

                eyvec = eyvec + Ps[is,js]*ymat0[:,:,js]
                epvec = epvec + Ps[is,js]*pmat0[:,:,js]
                evvec = evvec + Ps[is,js]*vmat0[:,:,js]

            end

            eyintp = LinearInterpolation((knotsm,knotsn),eyvec,extrapolation_bc=Line())
            epintp = LinearInterpolation((knotsm,knotsn),epvec,extrapolation_bc=Line())
            evintp = LinearInterpolation((knotsm,knotsn),evvec,extrapolation_bc=Line())

            @inbounds for im in 1:nm

                mpast = knotsm[im]

                @inbounds for in in 1:nn

                    npast = knotsn[in]

                    res = nlsolve(x->foc2(x,s0,mpast,npast,eyintp,epintp,bet,kap,invlam,sig,rstar),[mpast])
                    m0 = res.zero[1]
                    n0 = 0.0

                    ey0 = eyintp(m0,n0)
                    ep0 = epintp(m0,n0)
                    y0 = (1.0/bet)*npast - kap*m0
                    p0 = kap*y0 + bet*ep0 # if kap=0 and ep=0, p0=0
                    i0 = sig*(ey0-y0) + ep0 + s0 # x1 = ex - (1/sig)*(r1-ep-snow)

                    if (i0<=-rstar)

                        res = nlsolve(x->focr2(x,s0,mpast,npast,eyintp,epintp,bet,kap,invlam,sig,rstar),[mpast,npast])
                        m0 = res.zero[1]
                        n0 = res.zero[2]

                        ep0 = epintp(m0,n0)
                        y0 = -n0 + (1.0/bet)*npast - kap*m0
                        p0 = kap*y0 + bet*ep0 # if kap=0 and ep=0, p0=0
                        i0 = -rstar #0.0

                    end

                    v0 = -y0^2 - invlam*p0^2 + bet*evintp(m0,n0)

                    ymat1[im,in,is] = y0
                    pmat1[im,in,is] = p0
                    imat1[im,in,is] = i0
                    mmat1[im,in,is] = m0
                    nmat1[im,in,is] = n0
                    vmat1[im,in,is] = v0

                end

            end

        end

        ydiff = maximum(abs.(ymat1-ymat0))
        pdiff = maximum(abs.(pmat1-pmat0))
        vdiff = maximum(abs.(vmat1-vmat0))
        diff = maximum([ydiff pdiff vdiff])
        iter = iter + 1
        println([iter ydiff pdiff vdiff])

        ymat0 = damp*ymat1 + (1.0-damp)*ymat0
        pmat0 = damp*pmat1 + (1.0-damp)*pmat0
        imat0 = damp*imat1 + (1.0-damp)*imat0
        mmat0 = damp*mmat1 + (1.0-damp)*mmat0
        nmat0 = damp*nmat1 + (1.0-damp)*nmat0
        vmat0 = damp*vmat1 + (1.0-damp)*vmat0

    end

end

function calcmpe(Gs,Ps,bet,kap,invlam,sig,rstar)

    sh = Gs[1]
    sl = Gs[2]
    ph = Ps[1,2]
    pl = Ps[2,2]

    Ax = zeros(Float64,6,6)
    Ax[1,:] = [sig*ph -(1-ph) 1 -sig*ph -ph 0]
    Ax[2,:] = [-kap 1-bet*(1-ph) 0 0 -bet*ph 0]
    Ax[3,:] = [-1 -kap*invlam 0 0 0 0]
    Ax[4,:] = [-sig*(1-pl) -(1-pl) 0 sig*(1-pl) -pl 0]
    Ax[5,:] = [0 -bet*(1-pl) 0 -kap 1-bet*pl 0]
    Ax[6,:] = [0 0 0 -1 -kap*invlam -1]
    bx = zeros(Float64,6,1)
    bx[1] = sh #-rstar
    bx[4] = rstar+sl
    xx = inv(Ax)*bx
    yh = xx[1]
    yl = xx[4]
    pih = xx[2]
    pil = xx[5]
    rh = xx[3]
    rl = -rstar

    Av = zeros(Float64,2,2)
    Av[1,:] = [1-bet*(1-ph) -bet*ph]
    Av[2,:] = [-bet*(1-pl) 1-bet*pl]
    bv = zeros(Float64,2,1)
    bv[1] = -yh^2 - invlam*pih^2
    bv[2] = -yl^2 - invlam*pil^2
    vv = inv(Av)*bv;
    vh = vv[1]
    vl = vv[2]
    xvec0 = [yh yl]'
    pvec0 = [pih pil]'
    rvec0 = [rh rl]'
    vvec0 = [vh vl]'

    return xvec0, pvec0, rvec0, vvec0

end

function foc2(x,s0,mpast,npast,eyintp,epintp,bet,kap,invlam,sig,rstar)

    m0 = x[1]
    n0 = 0.0

    ey0 = eyintp(m0,n0)
    ep0 = epintp(m0,n0)

    y0 = (1.0/bet)*npast - kap*m0
    p0 = kap*y0 + bet*ep0 # if kap=0 and ep=0, p0=0

    f = -invlam*p0 + (1.0/bet/sig)*npast + m0 - mpast # if invlam=0, m0=mpast-(1/bet/sig)*npast

    return f

end

function focr2(x,s0,mpast,npast,eyintp,epintp,bet,kap,invlam,sig,rstar)

    m0 = x[1]
    n0 = x[2]

    ey0 = eyintp(m0,n0)
    ep0 = epintp(m0,n0)

    y0 = -n0 + (1.0/bet)*npast - kap*m0
    p0 = kap*y0 + bet*ep0 # if kap=0 and ep=0, p0=0
    i0 = -rstar

    f = zeros(2)
    f[1] = -invlam*p0 + (1.0/bet/sig)*npast + m0 - mpast # if invlam=0, m0=mpast-(1/bet/sig)*npast
    f[2] = -y0 + ey0 - (1.0/sig)*(i0-ep0-s0)

    return f

end

main()
