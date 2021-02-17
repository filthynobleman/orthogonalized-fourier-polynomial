function CTilde = eigprod_funmap(C, Phi0, Order)
%CTilde = EIGPROD_FUNMAP(C, Phi0) Computes the extension of the functional
%map to eigenproducts. Here, C is the functional map and Phi0 is the value
%of the constant function on the source manifold.
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 24 September 2020

    if length(Phi0) > 1
        Phi0 = Phi0(1);
    end
    
    if nargin < 3
        N = 2;
    else
        N = Order;
    end
    
    K1 = size(C, 1) - 1;
    K2 = size(C, 2) - 1;
    
    
    % Compute CTilde
    Ctilde11 = C;
    Ctilde12 = Phi0 .* kron(C(1, 2:end), C(1:end, 2:end)) + ...
               Phi0 .* kron([zeros(1, K2); C(2:end, 2:end)], C(1, 2:end));
    Ctilde21 = zeros(K1 * K1, K2 + 1);
    Ctilde22 = kron(C(2:end, 2:end), C(2:end, 2:end));

    CTilde = [Ctilde11, Ctilde12;
              Ctilde21, Ctilde22];
          
    % Remove repeated indices
%     Remove1 = setdiff(1:size(CTilde, 1), sub_index(2, K1));
    Remove1 = [];
    for i = 1:K1
        for j = i+1:K1
            I1 = (i - 1) * K1 + j + K1 + 1;
            I2 = (j - 1) * K1 + i + K1 + 1;
            CTilde(I1, :) = CTilde(I1, :) + CTilde(I2, :);
            Remove1(end + 1) = I2;
        end
    end
    Remove2 = setdiff(1:size(CTilde, 2), sub_index(2, K2));
    CTilde(Remove1, :) = [];
    CTilde(:, Remove2) = [];
    
    % The code above is equivalend to the folloqing piece of code, but it
    % is way more clean. Also, it is easily generalized to higher orders.
%           
%     remove = [];
%     for i = 1:K1
%         for j = i+1:K1
%             remove(end + 1) = (j - 1) * K1 + i;
%         end
%     end
%     CTilde(remove + K1 + 1, :) = [];
%     remove = [];
%     for i = 1:K2
%         for j = i+1:K2
%             remove(end + 1) = (j - 1) * K2 + i;
%         end
%     end
%     CTilde(:, remove + K2 + 1) = [];

end

