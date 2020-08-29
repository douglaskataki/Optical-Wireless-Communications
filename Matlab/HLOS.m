%% Função definida:
% xl,yl,zl -> dimensões da sala
% Area,FOV -> propriedades do Led Rx
% n(ordem do lambertiano)-> propriedade do Led Tx
% Nx,Ny -> propriedades do grid para a sala
% N_Tx, N_Rx -> normal dos planos para Tx e Rx
% Pos_Tx -> Posição do Tx
% ---------------------------------
% H = HLOS(Area,FOV,n,N,Nx,Ny,x,y,N_Tx,Pos_Tx,N_Rx)

function H = HLOS(Area,FOV,n,N,Nx,Ny,x,y,N_Tx,Pos_Tx,N_Rx)
    h = zeros(N);  
    for ii=1:Nx
        for jj=1:Ny
            % plano receptor
            Pos_Rx = [x(ii) y(jj) 0];
            % vetor distancia entre Tx e Rx
            vTxRx = Pos_Tx-Pos_Rx;
            d = sqrt(sum((vTxRx).^2));
            % ângulo entre Tx e Rx
            cosphi = abs(dot(vTxRx,N_Rx))/d;
            phi_los = acosd(cosphi);
            if(phi_los<=FOV)
                costheta = abs(dot(vTxRx,N_Tx))/d;
                h(ii,jj) = h (ii,jj)+(n+1)*Area*costheta^n/(2*pi*d^2);
            end
        end
    end
    H = h;
end