%%% Calculates the analytical solution


Uean = cell(Nx,Nz,ceil(tmax/dt));
Uhan = cell(Nx,Nz,ceil(tmax/dt));

for s=1:tmax/dt
    t=(s-1)*dt; 
 
    for i=1:Nx
        xc = dx/2 + (i-1)*dx;
            for k=1:Nz
                zc = dz/2 + (k-1)*dz;
                EH = eh_te10(xc,zc,t,a,A10,om,mu,eps,beta);
    
                Uean{i,k,s} = EH(1,1);
                Uhan{i,k,s} = EH(2:3,1);
            end
    end
end

