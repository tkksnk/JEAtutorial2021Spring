# using NLsolve
# using Interpolations
include("rw2.jl")

function main()

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

    m = Model(bet,kap,phi,rstar,zeta,sH,sL,pH,pL)
    # policy function iteration
    Ns = 2
    N = 101
    mmax = 0.0
    mmin = -6.0

    ti(m,mmin,mmax,Ns,N,damp)
    # ti(bet,kap,phi,rstar,zeta,pH,pL,sH,sL,mmin,mmax,Ns,N,damp)

end

main()
