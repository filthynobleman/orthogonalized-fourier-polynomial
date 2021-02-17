function CTilde = Cp_of_C(C, order, k1, k2, Phi0)
%CTilde = CP_OF_C(C, order, k1, k2, Phi0) Compute the functional map for
%eigenproducts of order N = 1, 2, 3. The matrix C maps the first k2
%non-constant eigenfunctions from the source shape onto the space of the
%first k1 non-constant eigenfunctions on the target shape. Phi0 is the
%constant value assumed by the constant eigenfunction on the source shape.
%
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 20 October 2020

    cN = Phi0;

    % building Cp, Bp
    CTilde=C; C0=CTilde(2:end, 2:end);
    sumpow1=1; sumpow2=1;
    for p=2:order
        % sumpow = 1 + k1 + k1² + ... + k1^(p-1)
        new_sumpow1=k1*sumpow1 +1;
        new_sumpow2=k2*sumpow2 +1;
        CTilde11 = CTilde;
        CTilde12 = cN * kron(CTilde(:, sumpow2+1:end), C(1,2:end)) + ...
                   [zeros(1, k2^p); kron([cN * CTilde(1, sumpow2+1:end); CTilde(2:sumpow1, sumpow2+1:end)], C0)];
        CTilde21 = zeros(k1^p, new_sumpow2);
        CTilde22 = kron(CTilde(sumpow1+1:end, sumpow2+1:end),C0);
        CTilde = [CTilde11, CTilde12;
                  CTilde21, CTilde22];
        sumpow1=new_sumpow1;
        sumpow2=new_sumpow2;
    end

    % Remove repeated indices
%     Remove1 = setdiff(1:size(CTilde, 1), sub_index(2, K1));
    if order > 1
        Remove1 = [];
        K1 = k1;
        K2 = k2;
        for i = 1:K1
            for j = i+1:K1
                I1 = (i - 1) * K1 + j + K1 + 1;
                I2 = (j - 1) * K1 + i + K1 + 1;
                CTilde(I1, :) = CTilde(I1, :) + CTilde(I2, :);
                Remove1(end + 1) = I2;
            end
        end
    end
    if order > 2
        for i = 1:K1
            for j = i:K1
                for k = j:K1
                    I1 = (i - 1) * K1^2 + (j - 1) * K1 + k + K1^2 + K1 + 1;    % i, j, k
                    I2 = (i - 1) * K1^2 + (k - 1) * K1 + j + K1^2 + K1 + 1;    % i, k, j
                    I3 = (j - 1) * K1^2 + (i - 1) * K1 + k + K1^2 + K1 + 1;    % j, i, k
                    I4 = (j - 1) * K1^2 + (k - 1) * K1 + i + K1^2 + K1 + 1;    % j, k, i
                    I5 = (k - 1) * K1^2 + (i - 1) * K1 + j + K1^2 + K1 + 1;    % k, i, j
                    I6 = (k - 1) * K1^2 + (j - 1) * K1 + i + K1^2 + K1 + 1;    % k, j, i

                    if i == j && j == k
                        % i = j = j only one entry, continue
                        assert(I1 == I2 && I2 == I3 && I3 == I4 && I4 == I5 && I5 == I6);
                        continue;
                    elseif i == j
                        % i = j ==> I1 = I3, I2 = I4, I5 = I6
                        assert(I1 == I3 && I2 == I4 && I5 == I6);
                        CTilde(I1, :) = CTilde(I1, :) + CTilde(I2, :) + CTilde(I5, :);
                        Remove1(end + 1) = I2;
                        Remove1(end + 1) = I5;
                    elseif i == k
                        % i = k ==> I1 == I6, I2 == I5, I3 == I4
                        assert(I1 == I6 && I2 == I5 && I3 == I4)
                        CTilde(I1, :) = CTilde(I1, :) + CTilde(I2, :) + CTilde(I3, :);
                        Remove1(end + 1) = I2;
                        Remove1(end + 1) = I3;
                    elseif j == k
                        % j = k ==> I1 == I2, I3 == I5, I4 == I6
                        assert(I1 == I2 && I3 == I5 && I4 == I6)
                        CTilde(I1, :) = CTilde(I1, :) + CTilde(I3, :) + CTilde(I4, :);
                        Remove1(end + 1) = I3;
                        Remove1(end + 1) = I4;
                    else
                        assert(length(unique([I1, I2, I3, I4, I5, I6])) == 6);
                        CTilde(I1, :) = CTilde(I1, :) + CTilde(I2, :) + CTilde(I3, :) + ...
                                        CTilde(I4, :) + CTilde(I5, :) + CTilde(I6, :);
                        Remove1(end + 1) = I2;
                        Remove1(end + 1) = I3;
                        Remove1(end + 1) = I4;
                        Remove1(end + 1) = I5;
                        Remove1(end + 1) = I6;
                    end
                end
            end
        end
    end
    Remove2 = setdiff(1:size(CTilde, 2), sub_index(order, K2));
    CTilde(Remove1, :) = [];
    CTilde(:, Remove2) = [];
