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

figure(1)
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

 
