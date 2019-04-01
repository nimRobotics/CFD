function [outputArg1,outputArg2] = bcsApply(wMat,pMat)
% applying bc on omega
    for j in range(ny+1):
        wMat[j,0] = 2*(pMat[0,j]-pMat[1,j])/(dx[0]**2)    # left wall
        wMat[j,nx] = 2*(pMat[nx,j]-pMat[nx-1,j])/(dx[0]**2) # right wall
    for i in range(nx+1):
        wMat[0,i] = 2*(pMat[i,0]-pMat[i,1])/(dy[0]**2)    # bottom wall
        wMat[ny,i] = 2*(pMat[i,ny]-pMat[i-1,ny])/(dy[0]**2)-(2*U)/dy[0] # top wall

    %applying BC, psi is zero at all the boundaries
    pMat[0,:]=0    # left wall
    pMat[nx,:]=0   # right wall
    pMat[:,0]=0    # bottom wall
    pMat[:,ny]=0   # top wall

    outputArg1 = flipud(wMat)
    outputArg2 = flipud(pMat))
end