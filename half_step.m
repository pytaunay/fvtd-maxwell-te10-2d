%%% Calculate the flux at n + 1/2 
for i=1:Nx
    for k=1:Nz
        Un = zeros(3,1);

        % Get my neighbors and cell surface
        NN_idx = cell(Ns,1);
        NN_idx{1} = [i-1,k,dz]'; 
        NN_idx{2} = [i+1,k,dz]'; 
        NN_idx{3} = [i,k-1,dx]'; 
        NN_idx{4} = [i,k+1,dx]'; 

        % Get my E-H field
        Ei = Ue{i,k}; 
        Hi = Uh{i,k}; 

        %%% Get the flux for each neighbor
        % Initialize flux vector
        F = zeros(3,1);
        for l=1:Ns
            idx = NN_idx{l};
            coord = idx(1:2);
            S = idx(3);
            
            %%% Boundary condition: perfect electric conductor
            if coord(1) <= 0 || coord(1) > Nx %|| coord(2) > Nz
                F = F + (Tpec*Apmat{l}*[Ei;Hi])*S;

            %%% Boundary condition: TE10 mode excitation
            elseif coord(2)<= 0 %|| coord(2) > Nz
                % Center of surface coordinate
                xs = dx*(i-1) + dx/2;
                if( coord(2) <= 0)
                    zs = 0;
                else
                    zs = L;
                end
                    
                EHs = eh_te10(xs,zs,t,a,A10,om,mu,eps,beta);
                F = F + (Apmat{l}*[Ei;Hi]+Ammat{l}*EHs)*S;

                %EHs = eh_te10(xs,zs,0,a,A10,om,mu,eps,beta);
                % No outgoing flux
                %F = F + (Ammat{l}*EHs)*S;

            %%% Absorbing BC -- no incoming flux
            elseif coord(2) > Nz
               F = F+ (Apmat{l}*[Ei;Hi])*S;
                El = Ue{coord(1),1};
                Hl = Uh{coord(1),1};
               %F = F+ (Apmat{l}*[Ei;Hi]+Ammat{l}*[El;Hl])*S;

            %%% General case
            else
                % Get the neighbor's fields
                El = Ue{coord(1),coord(2)};
                Hl = Uh{coord(1),coord(2)};
                
                if(coord(1) > i || coord(2) > k)
                    F = F + (Apmat{l}*[Ei;Hi] + Ammat{l}*[El;Hl])*S;
                elseif(coord(1) < i || coord(2) < k)
                    F = F + (Apmat{l}*[El;Hl] + Ammat{l}*[Ei;Hi])*S;
                end

                % Calculate the flux
                %F = F + (Apmat{l}*[Ei;Hi] + Ammat{l}*[El;Hl])*S;

            end
        end % for l=1:Ns

        % Take a step
        Un = cat(1,Ei,Hi) - 1/dS*dt*F; 
        Uen{i,k} = Un(1,1); 
        Uhn{i,k} = Un(2:3,1); 
    end % for k=1:Nz
end % for i=1:Nx 
