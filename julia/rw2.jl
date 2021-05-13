# Reifschneider-Williams rule
using NLsolve
using Interpolations

# struct Model{TI<:Integer, TF<:AbstractFloat}
#     rstar::TF   # pH=0のときの、定常状態での名目金利の値
#     bet::TF     # 割引率
#     phi::TF     # テイラー係数
#     kap::TF     # フィリップス曲線の傾き
#     zeta::TF     # フィリップス曲線の傾き
#     sH::TF      # 状態Hでの自然利子率の値
#     sL::TF      # 状態Lでの自然利子率の値
#     pH::TF      # 危機が起こる確率
#     pL::TF      # 危機の継続確率
#     # maxiter::TI # 繰り返し回数の最大値
#     # tol::TF     # 許容誤差
# end

# function main()
#
#     damp = 1.0
#
#     bet = 0.9925
#     kap = 0.01
#     phi = 5.0
#     rstar = (1.0/bet-1.0)*100
#     zeta = 1.0
#
#     pH = 0.0
#     pL = 0.8
#     sH = rstar
#     sL = rstar-2.5
#
#     m = Model(bet,kap,phi,rstar,zeta,sH,sL,pH,pL)
#     # policy function iteration
#     Ns = 2
#     N = 101
#     mmax = 0.0
#     mmin = -6.0
#
#     ti(m,mmin,mmax,Ns,N,damp)
#     # ti(bet,kap,phi,rstar,zeta,pH,pL,sH,sL,mmin,mmax,Ns,N,damp)
#
# end

function ti(m,mmin,mmax,Ns,N,damp)

    bet = m.bet
    kap = m.kap
    phi = m.phi
    rstar = m.rstar
    zeta = m.zeta
    pH = m.pH
    pL = m.pL
    sH = m.sH
    sL = m.sL

    Gs = [sH; sL]
    Ps = [1.0-pH pH; 1.0-pL pL]

    # mgrid = collect(LinRange(mmin,mmax,N))
    mgrid = range(mmin,stop=mmax,length=N)

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

    ymat0 = ones(Float64,N,1)*[yH yL]
    pmat0 = ones(Float64,N,1)*[piH piL]
    imat0 = ones(Float64,N,1)*[iH 0]
    # ymat0 = zeros(Float64,N,Ns)
    # pmat0 = zeros(Float64,N,Ns)
    # imat0 = zeros(Float64,N,Ns)
    mmat0 = zeros(Float64,N,Ns)
    nmat0 = zeros(Float64,N,Ns) # nortional rate
    ymat1 = copy(ymat0)
    pmat1 = copy(pmat0)
    imat1 = copy(imat0)
    mmat1 = copy(mmat0)
    nmat1 = copy(nmat0)

    # crit = 1e-4
    diff = 1e+4
    iter = 0

    while (diff>tol && iter<maxiter)

        @inbounds for is in 1:Ns

            s0 = Gs[is]

            eyvec = zeros(Float64,N)
            epvec = zeros(Float64,N)

            @inbounds for js in 1:Ns

                eyvec = eyvec + Ps[is,js]*ymat0[:,js]
                epvec = epvec + Ps[is,js]*pmat0[:,js]

            end

            eyintp = LinearInterpolation(mgrid,eyvec,extrapolation_bc=Line())
            epintp = LinearInterpolation(mgrid,epvec,extrapolation_bc=Line())

            @inbounds for im in 1:N

                mpast = mgrid[im]

                # solve for m0
                res = nlsolve(x->rw_eqm2(x,mpast,epintp,rstar,phi,zeta),[mpast])
                m0 = res.zero[1]

                ey0 = eyintp(m0)
                ep0 = epintp(m0)
                is0 = rstar + phi*ep0
                i0 = max(is0+mpast,0.0)
                y0 = ey0 - (i0 - ep0 - s0)
                p0 = kap*y0 + bet*ep0
                if (is0>0.0)
                    m0 = mpast - (i0-is0)
                else
                    m0 = mpast - zeta*(0-is0)
                end

                ymat1[im,is] = y0
                pmat1[im,is] = p0
                imat1[im,is] = i0
                mmat1[im,is] = m0
                nmat1[im,is] = is0

            end

        end

        ydiff = maximum(abs.(ymat1-ymat0))
        pdiff = maximum(abs.(pmat1-pmat0))
        diff = maximum([ydiff pdiff])
        iter = iter + 1
        println([iter diff])

        ymat0 = damp*ymat1 + (1.0-damp)*ymat0
        pmat0 = damp*pmat1 + (1.0-damp)*pmat0
        imat0 = damp*imat1 + (1.0-damp)*imat0
        mmat0 = damp*mmat1 + (1.0-damp)*mmat0
        nmat0 = damp*nmat1 + (1.0-damp)*nmat0

    end

    return ymat0, pmat0, imat0, mmat0, nmat0

end

function rw_eqm2(x,mpast,epintp,rstar,phi,zeta)

    m0 = x[1]
    epi = epintp(m0)
    istar = rstar + phi*epi
    i0 = max(istar+mpast,0.0)

    if (istar<0.0)
        f = -m0 + mpast - zeta*(0.0-istar)
    else
        f = -m0 + mpast - (i0-istar)
    end

    return f

end

# main()
