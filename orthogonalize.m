function [Ortho, Coeff, SubIdx] = orthogonalize(Area, Basis, Tol, Adjust)
%Ortho = ORTHOGONALIZE(Area, Basis) Applies the Gram-Schmidt algorithm to
%orthogonalize the given basis with an inner product depending on the given
%areas. If the basis vectors are not linearly independent, the algorithms
%discards each linearly dependent vector.
%
%Ortho = ORTHOGONALIZE(Area, Basis, Tol) Set the tolerance for the zero
%vector detection. Default is 1e-9;
%
%Ortho = ORTHOGONALIZE(Area, Basis, Tol, Adjust) Determines if the
%algorithm should try to adjust the vectors when they numerically turn out
%to be linearly dependent from the previous ones. If true, each iteration
%is applied to a vector multiple times, until it eventually becomes
%othogonal. Default is true.
%
%[Ortho, Coeff] = ORTHOGONALIZE(--) The algorithm also returns a matrix of
%coefficients such that Ortho * Coeff = Basis. This is particularly useful
%when Adjust = true, since the coefficients are not updated in the
%sdjustment and provides way more stable results.
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 23 September 2020

    if nargin < 4
        Adjust = true;
    end
    if nargin < 3
        Tol = 1e-9;
    end

    Ortho = Basis(:, 1) ./ sqrt(Basis(:, 1)' * Area * Basis(:, 1));
    if nargout > 1
        Coeff = sqrt(Basis(:, 1)' * Area * Basis(:, 1));
    end
    if nargout > 2
        SubIdx = 1;
    end
    
    for i = 2:size(Basis, 2)
        NewProj = Ortho' * (Area * Basis(:, i));
        New = Basis(:, i) - Ortho * NewProj;
        Norm = sqrt(New' * Area * New);
        if Norm < Tol
            if nargout > 1
                Coeff = [Coeff, NewProj];
            end
            continue;
        end
        New = New ./ Norm;
        if Adjust
            while any((Ortho'* (Area * New)) > Tol)
                CurProj = Ortho'* (Area * New);
                New = New - Ortho * CurProj;
%                 % Approach Wrong (Debug)
%                 Norm = sqrt(New'* Area * New);
%                 New = New ./ Norm;
                % Approach 1 (Release)
                CurNorm = sqrt(New'* Area * New);
                New = New ./ CurNorm;
%                 % Approach 2 (Discarded)
%                 CurNorm = sqrt(New'* Area * New);
%                 New = New ./ CurNorm;
%                 Norm = Norm * CurNorm;
%                 % Approach 3 (Discarded)
%                 CurNorm = sqrt(New'* Area * New);
%                 New = New ./ CurNorm;
%                 NewProj = NewProj + CurProj * Norm;
%                 Norm = Norm * CurNorm;
            end
        end
        Ortho = [Ortho, New];
        NewProj(abs(NewProj) < Tol) = 0;
        if nargout > 1
            Coeff = [Coeff, NewProj; 
                     zeros(1, size(Coeff, 2)), Norm];
        end
        if nargout > 2
            SubIdx(end + 1) = i;
        end
    end

end

