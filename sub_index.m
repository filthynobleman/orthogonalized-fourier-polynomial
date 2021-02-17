function [subs_order] = sub_index(order, k)


    PolyPhi = 0:k;
    Basis = 1:k;
    Prods = Basis;
    for i = 2:order
        tmp1 = repmat(Basis, 1, k^(i - 1));
        tmp2 = reshape(repmat(Prods, k, 1), i - 1, k^i);
        Prods = [tmp1; tmp2];
        PolyPhi = [PolyPhi; zeros(1, length(PolyPhi))];
        PolyPhi = [PolyPhi, Prods];
    end
    PolyPhi = sort(PolyPhi, 1);
    [PolyPhi, subs_order, ~] = unique(PolyPhi', 'rows');
    

%     a=k+order-1; b=order;
%     choose=zeros(a, b +1);
%     choose(:, 1)=ones(a, 1); choose(1, 1 +1)=1;
%     for i=1:b
%         for j=2:a
%             choose(j, i +1)=choose(j-1, i +1)+ choose(j-1, i-1 +1);
%         end
%     end
% 
%     subs_order=1:(1+k);
% 
%     for ord=2:order
% 
%         poss=cell(k-1, ord);
%         for j=1:(k-1)
%             poss{j, ord+1 -1}=[];
%         end
%         for i=ord:-1:2
%              poss{1, i -1}=[choose(k+i-1-1, i-1 +1)];
%             for j=2:(k-1)
%         %        poss{j, i -1}=choose(k+i-j-1, i-1 +1)+[poss{j-1, i+1 -1}; poss{j-1, i -1}];
%                 poss{j, i -1}=[-choose(k+i-j-2, i-1 +1)+poss{j, i+1 -1}; choose(k+i-j-1, i-1 +1)+poss{j-1, i -1}];
%             end
%         end
% 
%         sumpowk=ones(ord,1);
%         for i=1:(ord-1)
%             sumpowk(i+1)=sumpowk(i)*k;
%         end
%         sumpowk=cumsum(sumpowk);
% 
%         vals=zeros(k-1, ord-1);
%         for i=2:ord
%             for j=1:(k-1)
%                 vals(j, i -1)=j*sumpowk(i-1);
%             end
%         end
% 
%         skip_pos=[]; skip_val=[];
%         for i=2:ord
%             for j=1:(k-1)
%                 skip_val=[skip_val; repmat(vals(j, i -1), length(poss{j, i -1}), 1)];
%                 skip_pos=[skip_pos; poss{j, i -1}];
%             end
%         end
% 
%         subs_order=[subs_order, subs_order(end)+cumsum(ones(choose(k+ord-1, ord +1), 1)+ accumarray(1+skip_pos, skip_val))'];
% 
%     end
end