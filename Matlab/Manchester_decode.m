%% Função definida:
% bits = os bits utilizados para o vetor codificado

function MD = Manchester_decode(manchester)
    N = length(manchester)/2;
    new_code = zeros(1,N);
    for ii=1:2:2*N-1
       if isequal([manchester(ii) manchester(ii+1)],[1 0])
           new_code((ii+1)/2)=0;
       else
           new_code((ii+1)/2)=1;
       end
    end
    MD = new_code;
end