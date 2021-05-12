using NLsolve
using Interpolations

function main()
    # Reifschneider-Williams rule
    damp = 1.0

    bet = 0.9925
    kap = 0.01
    phi = 5.0
    rstar = (1/bet-1)*100
    zeta = 1.0

    pH = 0.0
    pL = 0.8
    sH = rstar
    sL = rstar-2.5

    # initial guess for zeta = 0
    # (semi-)analytical
    A = [-1.0+(1.0-pH) pH -(phi-1.0)*(1.0-pH) -(phi-1.0)*pH;
        kap 0.0 -1.0+bet*(1.0-pH) bet*pH;
        (1.0-pL) -1.0+pL (1.0-pL) pL;
        0.0 kap bet*(1.0-pL) -1.0+bet*pL]
    b = [rstar-sH;0.0;-sL;0.0]
    x = A\b

    yH  = x[1]
    yL  = x[2]
    piH = x[3]
    piL = x[4]
    # check>0
    iH = rstar + phi*((1-pH)*piH + pH*piL)
    println(iH)

    N = 101
    Ns = 2
    mmax = 0.0
    mmin = -6.0
    # mgrid = collect(LinRange(mmin,mmax,N))
    mgrid = range(mmin,stop=mmax,length=N)

    Gs = [sH; sL]
    Ps = [1.0-pH pH; 1.0-pL pL]

    invlam = 1/0.003
    sig = 1.0
    xvec0, pvec0, rvec0, vvec0 = calcmpe(Gs,Ps,bet,kap,invlam,sig,rstar)
    println([xvec0[2] pvec0[2]])

    # policy function iteration
    ymat0 = ones(N,1)*[yH yL]
    pmat0 = ones(N,1)*[piH piL]
    imat0 = ones(N,1)*[iH 0]
    ymat0 = zeros(Float64,N,Ns)
    pmat0 = zeros(Float64,N,Ns)
    imat0 = zeros(Float64,N,Ns)
    mmat0 = zeros(Float64,N,Ns)
    nmat0 = zeros(Float64,N,Ns) # nortional rate
    ymat1 = copy(ymat0)
    pmat1 = copy(pmat0)
    imat1 = copy(imat0)
    mmat1 = copy(mmat0)
    nmat1 = copy(nmat0)

    # crit = 1e-4
    # diff = 1e+4
    # iter = 0
    #
    # while (diff>crit)
    #
    #     @inbounds for is in 1:Ns
    #
    #         s0 = Gs[is]
    #
    #         eyvec = zeros(Float64,N)
    #         epvec = zeros(Float64,N)
    #
    #         @inbounds for js in 1:Ns
    #
    #             eyvec = eyvec + Ps[is,js]*ymat0[:,js]
    #             epvec = epvec + Ps[is,js]*pmat0[:,js]
    #
    #         end
    #
    #         eyintp = LinearInterpolation(mgrid,eyvec,extrapolation_bc=Line())
    #         epintp = LinearInterpolation(mgrid,epvec,extrapolation_bc=Line())
    #
    #         @inbounds for im in 1:N
    #
    #             mpast = mgrid[im]
    #
    #             # solve for m0
    #             # res = nlsolve(x->rw_eqm(x,mpast,mgrid,epvec,rstar,phi,zeta),[mpast])
    #             res = nlsolve(x->rw_eqm2(x,mpast,epintp,rstar,phi,zeta),[mpast])
    #             m0 = res.zero[1]
    #
    #             ey0 = eyintp(m0)
    #             ep0 = epintp(m0)
    #             # ey0 = intf1(mgrid,eyvec,m0)
    #             # ep0 = intf1(mgrid,epvec,m0)
    #             is0 = rstar + phi*ep0
    #             i0 = max(is0+mpast,0.0)
    #             y0 = ey0 - (i0 - ep0 - s0)
    #             p0 = kap*y0 + bet*ep0
    #             if (is0>0.0)
    #                 m0 = mpast - (i0-is0)
    #             else
    #                 m0 = mpast - zeta*(0-is0)
    #             end
    #
    #             ymat1[im,is] = y0
    #             pmat1[im,is] = p0
    #             imat1[im,is] = i0
    #             mmat1[im,is] = m0
    #             nmat1[im,is] = is0
    #
    #         end
    #
    #     end
    #
    #     ydiff = maximum(abs.(ymat1-ymat0))
    #     pdiff = maximum(abs.(pmat1-pmat0))
    #     # idiff = maximum(abs.(imat1-imat0))
    #     diff = maximum([ydiff pdiff])
    #     iter = iter + 1
    #     println([iter diff])
    #
    #     ymat0 = damp*ymat1 + (1.0-damp)*ymat0
    #     pmat0 = damp*pmat1 + (1.0-damp)*pmat0
    #     imat0 = damp*imat1 + (1.0-damp)*imat0
    #     mmat0 = damp*mmat1 + (1.0-damp)*mmat0
    #     nmat0 = damp*nmat1 + (1.0-damp)*nmat0
    #
    # end

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
    bx[1] = sh-rstar
    bx[4] = sl
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

main()
