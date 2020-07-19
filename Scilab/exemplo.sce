clc;
clear;
// Função definida:
// xl,yl,zl -> dimensões da sala
// Area,FOV -> propriedades do Led Rx
// n(ordem do lambertiano)-> propriedade do Led Tx
// Nx,Ny -> propriedades do grid para a sala
// N_Tx, N_Rx -> normal dos planos para Tx e Rx
// Pos_Tx -> Posição do Tx
function H = HLOS(Area,FOV,n,N,Nx,Ny,x,y,N_Tx,Pos_Tx,N_Rx)
    h = zeros(N,N);  
    for ii=1:Nx
        for jj=1:Ny
            // plano receptor
            Pos_Rx = [x(ii) y(jj) 0];
            // vetor distancia entre Tx e Rx
            vTxRx = Pos_Tx-Pos_Rx;
            d = sqrt(sum((vTxRx).^2));
            // ângulo entre Tx e Rx
            cosphi = abs(sum(vTxRx.*N_Rx))/d;
            phi_los = acosd(cosphi);
            if(phi_los<=FOV)
                costheta = abs(sum(vTxRx.*N_Tx))/d;
                h(ii,jj) = h (ii,jj)+(n+1)*Area*costheta^n/(2*%pi*d^2);
            end
        end
    end
    H = h;
endfunction

// Função definida:
// xl,yl,zl,rho(coef reflexão) -> dimensões da sala
// n(ordem do lambertiano)-> propriedade do Led Tx
// Area,FOV-> propriedades do Led Rx
// Nx,Ny,Nz,x,y,z -> propriedades do grid para a sala (Quantidade e vetor)
// N_Tx, N_Rx -> normal dos planos para Tx e Rx
// Pos_Tx -> Posição do Tx
function HN = HNLOS(xl,yl,zl,rho,Area,FOV,n,Nx,Ny,Nz,x,y,z,N_Tx,Pos_Tx,N_Rx)
    N = max([Nx Ny Nz]);
    h1 = zeros(N,N);
    h2 = zeros(N,N);
    h3 = zeros(N,N);
    h4 = zeros(N,N);
    for ii=1:Nx
        for jj=1:Ny
            Pos_Rx = [x(ii) y(jj) 0];
            for kk=1:Ny
            //////////////// plano X=0 //////////////////
            // fator de área
            dA = zl*yl/(Ny*Nz);
            // normal da parede
            np = [1 0 0];
                for ll=1:Nz
                    // ponto na parede Z,Y - Wall Point
                    WP = [0 y(kk) z(ll)];
                    vTxWP = Pos_Tx - WP; 
                    // distancia do TX para a parede Z,Y(1)
                    D1 = sqrt(sum((vTxWP).^2));
                    // distancia TX e Rx
                    // angulos Tx e incidencia na parede
                    cos_phi = abs(sum(N_Tx.*vTxWP))/D1;
                    cos_alpha = abs(sum(vTxWP.*np))/D1;
                    // distancia do WP para o Rx
                    vWPRx = WP-Pos_Rx;
                    D2 = sqrt(sum((vWPRx).^2));
                    // angulos Rx e reflexão
                    cos_psi = abs(sum(vWPRx.*N_Rx))/D2;
                    cos_beta = abs(sum(vWPRx.*np))/D2;
                    if abs(acosd(cos_psi))<=FOV
                        h1(ii,jj) = h1(ii,jj)+(n+1)*Area*rho*dA*cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*%pi^2*D1^2*D2^2);
                    end
                end
            end
    
            //////////////// plano Y=0 //////////////////
            // fator de área
            dA = zl*xl/(Nx*Nz);
            // normal da parede
            np = [0 1 0];
            for kk=1:Nx
                for ll=1:Nz
                    // ponto na parede Z,Y - Wall Point
                    WP = [x(kk) 0 z(ll)];
                    vTxWP = Pos_Tx - WP; 
                    // distancia do TX para a parede Z,Y(1)
                    D1 = sqrt(sum((vTxWP).^2));
                    // distancia TX e Rx
                    // angulos Tx e incidencia na parede
                    cos_phi = abs(sum(N_Tx.*vTxWP))/D1;
                    cos_alpha = abs(sum(vTxWP.*np))/D1;
                    // distancia do WP para o Rx
                    vWPRx = WP-Pos_Rx;
                    D2 = sqrt(sum((vWPRx).^2));
                    // angulos Rx e reflexão
                    cos_psi = abs(sum(vWPRx.*N_Rx))/D2;
                    cos_beta = abs(sum(vWPRx.*np))/D2;
                    if abs(acosd(cos_psi))<=FOV
                        h2(ii,jj) = h2(ii,jj)+(n+1)*Area*rho*dA*cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*%pi^2*D1^2*D2^2);
                    end
                end
            end
            //////////////// plano X=xl //////////////////
            // fator de área
            dA = zl*yl/(Ny*Nz);
            // normal da parede
            np = [-1 0 0];
            for kk=1:Ny
                for ll=1:Nz
                    // ponto na parede Z,Y - Wall Point
                    WP = [xl y(kk) z(ll)];
                    vTxWP = Pos_Tx - WP; 
                    // distancia do TX para a parede Z,Y(1)
                    D1 = sqrt(sum((vTxWP).^2));
                    // distancia TX e Rx
                    // angulos Tx e incidencia na parede
                    cos_phi = abs(sum(N_Tx.*vTxWP))/D1;
                    cos_alpha = abs(sum(vTxWP.*np))/D1;
                    // distancia do WP para o Rx
                    vWPRx = WP-Pos_Rx;
                    D2 = sqrt(sum((vWPRx).^2));
                    // angulos Rx e reflexão
                    cos_psi = abs(sum(vWPRx.*N_Rx))/D2;
                    cos_beta = abs(sum(vWPRx.*np))/D2;
                    if abs(acosd(cos_psi))<=FOV
                        h3(ii,jj) = h3(ii,jj)+(n+1)*Area*rho*dA*cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*%pi^2*D1^2*D2^2);
                    end
                end
            end        
            //////////////// plano Y=yl //////////////////
            // fator de área
            dA = zl*xl/(Nx*Nz);
            // normal da parede
            np = [0 -1 0];
            for kk=1:Nx
                for ll=1:Nz
                    // ponto na parede Z,Y - Wall Point
                    WP = [x(kk) yl z(ll)];
                    vTxWP = Pos_Tx - WP; 
                    // distancia do TX para a parede Z,Y(1)
                    D1 = sqrt(sum((vTxWP).^2));
                    // distancia TX e Rx
                    // angulos Tx e incidencia na parede
                    cos_phi = abs(sum(N_Tx.*vTxWP))/D1;
                    cos_alpha = abs(sum(vTxWP.*np))/D1;
                    // distancia do WP para o Rx
                    vWPRx = WP-Pos_Rx;
                    D2 = sqrt(sum((vWPRx).^2));
                    // angulos Rx e reflexão
                    cos_psi = abs(sum(vWPRx.*N_Rx))/D2;
                    cos_beta = abs(sum(vWPRx.*np))/D2;
                    if abs(acosd(cos_psi))<=FOV
                        h4(ii,jj) = h4(ii,jj)+(n+1)*Area*rho*dA*cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*%pi^2*D1^2*D2^2);
                    end
                end
            end
        end
    end
    HN = h1+h2+h3+h4;
endfunction

// Propriedades no Detector
// Potencia LED
PLed = 20; // mW
// campo de visao do fotodetector
FOV = 60; // graus
// area do fotodetector
Area = 5.8E-6; // m^2
// semi-angulo do Tx para meia potencia
theta = 60; // graus
// n Ordem dem Emissao Lambertiana 
n = -log10(2)/log10(cosd(theta));
// Ganho optico no concentrador
index = 1.5;
G_Con = index^2/sin(FOV)^2;
TS = 1;
// Posicoes Tx e Rx
// tamanhos iniciais da sala
xl = 5; yl = 5; zl = 3; 
// altura do Rx
h = zl;
// Posicao Tx
x_Tx = xl/2; y_Tx = yl/2; z_Tx = zl;
Pos_Tx = [x_Tx y_Tx z_Tx];
// normaliza o locais para os Rx
M = 10;
Nx = xl*M; Ny = yl*M; Nz = round(h*M);
x = linspace(0,xl,Nx);
y = linspace(0,yl,Ny);    
z = linspace(0,zl,Nz);
N = max([Nx Ny Nz]);
// base para os Rx
N_Tx = [0 0 -1];
N_Rx = -N_Tx;

// Matrizes de Rotações
// ângulo de rotação, lembrar de colocar em graus
// angulosTx = [0 0 0];
// angulosRx = [0 0 0];

// Tx
// rotx = [1 0 0; 0 cosd(angulosTx(1)) -sind(angulosTx(1)); 0 sind(angulosTx(1)) cosd(angulosTx(1))];
// roty = [cosd(angulosTx(2)) 0 sind(angulosTx(2)); 0 1 0; -sind(angulosTx(2)) 0 cosd(angulosTx(2))];
// rotz = [cos(angulosTx(3)) -sin(angulosTx(3)) 0; sin(angulosTx(3)) cos(angulosTx(3)) 0; 0 0 1];
// N_Tx = rotz*roty*rotx*N_Tx;

// Rx
// rotx = [1 0 0; 0 cosd(angulosRx(1)) -sind(angulosRx(1)); 0 sind(angulosRx(1)) cosd(angulosRx(1))];
// roty = [cosd(angulosRx(2)) 0 sind(angulosRx(2)); 0 1 0; -sind(angulosRx(2)) 0 cosd(angulosRx(2))];
// rotz = [cos(angulosRx(3)) -sin(angulosRx(3)) 0; sin(angulosRx(3)) cos(angulosRx(3)) 0; 0 0 1];
// N_Rx = rotz*roty*rotx*N_Rx;

// Calculo da Resposta do Canal
H = HLOS(Area,FOV,n,N,Nx,Ny,x,y,N_Tx,Pos_Tx,N_Rx);
// potencia de chegada
P_Rx = PLed*G_Con.*H;
// para dBm
P_Rx_dBm = 10*log10(P_Rx);

subplot(1,2,1)
plot3d(x,y,P_Rx_dBm);
/*
h = gce();
h.color_flag=1;
h.color_mode=-2
h.color_flag=2;
h.color_mode = -1
f=gcf();
f.color_map=hotcolormap(512);
c=[1:400,1:400];
TL.color = [c;c+1;c+2;c+3];
h.data = TL;
h.color_flag=3;
*/
xtitle("HLOS - Led");
xlabel("x(m)");
ylabel("y(m)");
zlabel("dBm");

rho = 0.8;
Hlos = HNLOS(xl,yl,zl,rho,Area,FOV,n,Nx,Ny,Nz,x,y,z,N_Tx,Pos_Tx,N_Rx);
P_Rx = PLed*G_Con.*(H+Hlos);
// para dBm
P_Rx_dBm = 10*log10(P_Rx);

subplot(1,2,2)
plot3d(x,y,P_Rx_dBm);
xtitle("HLOS+1stHNLOS - Led");
xlabel("x(m)");
ylabel("y(m)");
zlabel("dBm");
