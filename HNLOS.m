%% Função definida:
% xl,yl,zl,rho(coef reflexão) -> dimensões da sala
% Area,FOV,n(ordem do lambertiano) -> propriedades do Led
% Nx,Ny,Nz,x,y,z -> propriedades do grid para a sala
% N_Tx, N_Rx -> normal dos planos para Tx e Rx
% Pos_Tx -> Posição do Tx
function H = HNLOS(xl,yl,zl,rho,Area,FOV,n,Nx,Ny,Nz,N,x,y,z,N_Tx,Pos_Tx,N_Rx)
h1 = zeros(N);
h2 = zeros(N);
h3 = zeros(N);
h4 = zeros(N);
for ii=1:Nx
    for jj=1:Ny
        Pos_Rx = [x(ii) y(jj) 0];
        for kk=1:Ny
        %%%%%%%%%%%%%%% plano X=0 %%%%%%%%%%%%%%%%%%
        % fator de área
        dA = zl*yl/(Ny*Nz);
        % normal da parede
        n1 = [1 0 0];
            for ll=1:Nz
                % ponto na parede Z,Y - Wall Point
                WP = [0 y(kk) z(ll)];
                vTxWP = Pos_Tx - WP; 
                % distancia do TX para a parede Z,Y(1)
                D1 = sqrt(sum((vTxWP).^2));
                % distancia TX e Rx
                % angulos Tx e incidencia na parede
                cos_phi = abs(dot(N_Tx,vTxWP))/D1;
                cos_alpha = abs(dot(vTxWP,n1))/D1;
                % distancia do WP para o Rx
                vWPRx = WP-Pos_Rx;
                D2 = sqrt(sum((vWPRx).^2));
                % angulos Rx e reflexão
                cos_psi = abs(dot(vWPRx,N_Rx))/D2;
                cos_beta = abs(dot(vWPRx,n1))/D2;
                if abs(acosd(cos_psi))<=FOV
                    h1(ii,jj) = h1(ii,jj)+(n+1)*Area*rho*dA*...
                    cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end

        %%%%%%%%%%%%%%% plano Y=0 %%%%%%%%%%%%%%%%%%
        % fator de área
        dA = zl*xl/(Nx*Nz);
        % normal da parede
        n2 = [0 1 0];
        for kk=1:Nx
            for ll=1:Nz
                % ponto na parede Z,Y - Wall Point
                WP = [x(kk) 0 z(ll)];
                vTxWP = Pos_Tx - WP; 
                % distancia do TX para a parede Z,Y(1)
                D1 = sqrt(sum((vTxWP).^2));
                % distancia TX e Rx
                % angulos Tx e incidencia na parede
                cos_phi = abs(dot(N_Tx,vTxWP))/D1;
                cos_alpha = abs(dot(vTxWP,n2))/D1;
                % distancia do WP para o Rx
                vWPRx = WP-Pos_Rx;
                D2 = sqrt(sum((vWPRx).^2));
                % angulos Rx e reflexão
                cos_psi = abs(dot(vWPRx,N_Rx))/D2;
                cos_beta = abs(dot(vWPRx,n2))/D2;
                if abs(acosd(cos_psi))<=FOV
                    h2(ii,jj) = h2(ii,jj)+(n+1)*Area*rho*dA*...
                    cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end
        %%%%%%%%%%%%%%% plano X=xl %%%%%%%%%%%%%%%%%%
        % fator de área
        dA = zl*yl/(Ny*Nz);
        % normal da parede
        n3 = [-1 0 0];
        for kk=1:Ny
            for ll=1:Nz
                % ponto na parede Z,Y - Wall Point
                WP = [xl y(kk) z(ll)];
                vTxWP = Pos_Tx - WP; 
                % distancia do TX para a parede Z,Y(1)
                D1 = sqrt(sum((vTxWP).^2));
                % distancia TX e Rx
                % angulos Tx e incidencia na parede
                cos_phi = abs(dot(N_Tx,vTxWP))/D1;
                cos_alpha = abs(dot(vTxWP,n3))/D1;
                % distancia do WP para o Rx
                vWPRx = WP-Pos_Rx;
                D2 = sqrt(sum((vWPRx).^2));
                % angulos Rx e reflexão
                cos_psi = abs(dot(vWPRx,N_Rx))/D2;
                cos_beta = abs(dot(vWPRx,n3))/D2;
                if abs(acosd(cos_psi))<=FOV
                    h3(ii,jj) = h3(ii,jj)+(n+1)*Area*rho*dA*...
                    cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end        
        %%%%%%%%%%%%%%% plano Y=yl %%%%%%%%%%%%%%%%%%
        % fator de área
        dA = zl*xl/(Nx*Nz);
        % normal da parede
        n4 = [0 -1 0];
        for kk=1:Nx
            for ll=1:Nz
                % ponto na parede Z,Y - Wall Point
                WP = [x(kk) yl z(ll)];
                vTxWP = Pos_Tx - WP; 
                % distancia do TX para a parede Z,Y(1)
                D1 = sqrt(sum((vTxWP).^2));
                % distancia TX e Rx
                % angulos Tx e incidencia na parede
                cos_phi = abs(dot(N_Tx,vTxWP))/D1;
                cos_alpha = abs(dot(vTxWP,n4))/D1;
                % distancia do WP para o Rx
                vWPRx = WP-Pos_Rx;
                D2 = sqrt(sum((vWPRx).^2));
                % angulos Rx e reflexão
                cos_psi = abs(dot(vWPRx,N_Rx))/D2;
                cos_beta = abs(dot(vWPRx,n4))/D2;
                if abs(acosd(cos_psi))<=FOV
                    h4(ii,jj) = h4(ii,jj)+(n+1)*Area*rho*dA*...
                    cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end
    end
end
H = h1+h2+h3+h4;
end