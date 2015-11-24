%%% Function flux_matrix
%%% Calculates the matrix from the split flux approach
%%% The matrix depends on the normal vector to a surface and the medium

function A = flux_matrix(n,mu,eps)

    nx = n(1,1);
   % ny = n(2,1);
    nz = n(3,1);

    c = 1/sqrt(mu*eps); 

    A = zeros(3,3);

    A = [c, -nz/eps, nx/eps;
        -nz/mu, c*nz^2, -c*nx*nz;
        nx/mu, -c*nx*nz, c*nx^2];

    %%% Result
    A = 1/2*A;

end
