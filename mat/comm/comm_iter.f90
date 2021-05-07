#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use mod_main
    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxCreateDoubleMatrix, mxGetPr
    mwPointer mxGetM, mxGetN
    mwPointer mxGetDimensions, mxGetNumberOfDimensions
    mwSize dims(3)
    mwPointer mxCreateNumericArray
    integer mxClassIDFromClassName

    ! [nmat1 mmat1 xmat1 pmat1 rmat1] = comm_iter2(Gs,Ps,knotsm,knotsn,mmat0,nmat0,xmat0,pmat0,bet,kap,invlam,sig,rstar)
    ! function [nmat1 xmat1 rmat1] = comm_iter(Gs,Ps,knotsn,nmat0,xmat0,bet,sig,rstar)
    real(8), allocatable :: Gs(:), Ps(:,:), knotsm(:), knotsn(:)
    real(8), allocatable :: mmat0(:,:,:), nmat0(:,:,:), xmat0(:,:,:), pmat0(:,:,:), rmat0(:,:,:), vmat0(:,:,:)
    real(8), allocatable :: mmat1(:,:,:), nmat1(:,:,:), xmat1(:,:,:), pmat1(:,:,:), rmat1(:,:,:), vmat1(:,:,:)
    integer ns, nm, nn, i
    real(8) bet, kap, invlam, sig, rstar

    integer, external :: mexPrintf
    integer k
    character(len=80) line

    call mxCopyPtrToInteger8(mxGetDimensions(prhs(5)), dims, mxGetNumberOfDimensions(prhs(5)))
    ns = dims(1)
    nm = dims(2)
    nn = dims(3)
    ! ns = mxGetM(prhs(1))
    ! nm = mxGetM(prhs(3))
    ! nn = mxGetM(prhs(4))
    allocate(Gs(ns),Ps(ns,ns),knotsm(nm),knotsn(nn))
    allocate(mmat0(ns,nm,nn),nmat0(ns,nm,nn),xmat0(ns,nm,nn),pmat0(ns,nm,nn),rmat0(ns,nm,nn),vmat0(ns,nm,nn))
    allocate(mmat1(ns,nm,nn),nmat1(ns,nm,nn),xmat1(ns,nm,nn),pmat1(ns,nm,nn),rmat1(ns,nm,nn),vmat1(ns,nm,nn))

    i = 1;
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Gs,ns); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ps,ns*ns); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),knotsm,nm); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),knotsn,nn); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),mmat0,ns*nm*nn); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),nmat0,ns*nm*nn); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),xmat0,ns*nm*nn); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),pmat0,ns*nm*nn); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),vmat0,ns*nm*nn); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),bet,1); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),kap,1); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),invlam,1); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),sig,1); i = i+1
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),rstar,1); i = i+1

    ! dims(1) = ns
    ! dims(2) = nm
    ! dims(3) = nn
    ! write(line,*) dims(1)
    ! k = mexPrintf(line//achar(10))
    ! write(line,*) dims(2)
    ! k = mexPrintf(line//achar(10))
    ! write(line,*) dims(3)
    ! k = mexPrintf(line//achar(10))
    call comm_iter(Gs,Ps,knotsm,knotsn,mmat0,nmat0,xmat0,pmat0,vmat0,bet,kap,invlam,sig,rstar,mmat1,nmat1,xmat1,pmat1,rmat1,vmat1)

    plhs(1) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(2) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(3) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(4) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(5) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(6) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)

    call mxCopyReal8ToPtr(mmat1,mxGetPr(plhs(1)),ns*nm*nn)
    call mxCopyReal8ToPtr(nmat1,mxGetPr(plhs(2)),ns*nm*nn)
    call mxCopyReal8ToPtr(xmat1,mxGetPr(plhs(3)),ns*nm*nn)
    call mxCopyReal8ToPtr(pmat1,mxGetPr(plhs(4)),ns*nm*nn)
    call mxCopyReal8ToPtr(rmat1,mxGetPr(plhs(5)),ns*nm*nn)
    call mxCopyReal8ToPtr(vmat1,mxGetPr(plhs(6)),ns*nm*nn)

    return

end subroutine mexFunction
