%% Função definida:
% bits = os bits utilizados para o vetor codificado

function MC = Manchester_code(bits)
    N = length(bits);
    new_code = zeros(1,2*N);
    for ii=1:2:2*N-1
       if bits((ii+1)/2)==1
           new_code(ii)=0;
           new_code(ii+1)=1;
       else
           new_code(ii)=1;
           new_code(ii+1)=0;
       end
    end
    MC = new_code;
end