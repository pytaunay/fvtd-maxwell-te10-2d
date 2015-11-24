[XG,ZG] = meshgrid(dx/2:dx:a,dz/2:dz:L);

figure();


Econtour = zeros(size(XG,1),size(XG,2));
Hcontour = zeros(size(XG,1),size(XG,2));

Econtour_a = zeros(size(XG,1),size(XG,2));
Hcontour_a = zeros(size(XG,1),size(XG,2));
for k=1:Nz
    for i=1:Nx
        Econtour(k,i) = sqrt(sum(abs(Ueall{i,k,1}).^2));
        Hcontour(k,i) = sqrt(sum(abs(Uhall{i,k,1}).^2));
        
        % Analytical
        Econtour_a(k,i) = sqrt(sum(abs(Uean{i,k,1}).^2));
        Hcontour_a(k,i) = sqrt(sum(abs(Uhan{i,k,1}).^2));
    end
end

subplot(2,2,1);
eCont=surf(XG,ZG,Econtour);
%axis([0,a,0,L,0,2.5e3]);
axis([0,a,0,L]);
subplot(2,2,2);
hCont=surf(XG,ZG,Hcontour);
axis([0,a,0,L,0,10]);

subplot(2,2,3);
eCont_a=surf(XG,ZG,Econtour_a);
axis([0,a,0,L]);
subplot(2,2,4);
hCont_a=surf(XG,ZG,Hcontour_a);
axis([0,a,0,L]);


pause;


for t=2:tmax/dt
    Econtour = zeros(size(XG,1),size(XG,2));
    Hcontour = zeros(size(XG,1),size(XG,2));
    
    Econtour_a = zeros(size(XG,1),size(XG,2));
    Hcontour_a = zeros(size(XG,1),size(XG,2));
    for k=1:Nz
        for i=1:Nx
        Econtour(k,i) = sqrt(sum(abs(Ueall{i,k,t}).^2));
        Hcontour(k,i) = sqrt(sum(abs(Uhall{i,k,t}).^2));
            
        Econtour_a(k,i) = sqrt(sum(abs(Uean{i,k,t}).^2));
        Hcontour_a(k,i) = sqrt(sum(abs(Uhan{i,k,t}).^2));
        end 
    end

    set(eCont,'ZData',Econtour);    
    set(hCont,'ZData',Hcontour);  
    
    set(eCont_a,'ZData',Econtour_a);
    set(hCont_a,'ZData',Hcontour_a);
pause(.1);
end
    
