function P2fine = cp_subpixelFine(P1,P2)
%CP_SUBPIXELFINE refine the coarse match coordinates in subpixel level
len = length(P1);
affmat = fitgeotrans(P1,P2,'projective');
P2pro = [P1 ones(len,1)] * affmat.T;
P2pro = P2pro(:,1:2) ./ P2pro(:,3);
devia_P = (P2 - P2pro).^2;
devia_P = sqrt(sum(devia_P,2));
max_Devia = max( devia_P );
iteration = 0;
P2fine = P2;
while max_Devia > 0.05 && iteration < 20
    iteration = iteration+1;
    fprintf('\nsubpixel iteration = %d\tmax_devia=%f\n',iteration,max_Devia);
    [~,index] = sort(devia_P);
    ind1 = round( 1/4 * length(index));
    P2fine(index(ind1:end),:)    =  P2pro(index(ind1:end),:);
    affmat = fitgeotrans(P1,P2fine,'projective');
    P2pro = [P1 ones(len,1)] * affmat.T;
    P2pro = P2pro(:,1:2) ./ P2pro(:,3);
    devia_P = (P2fine - P2pro).^2;
    devia_P = sqrt(sum(devia_P,2));
    max_Devia = max( devia_P );
end
fprintf('\nsubpixel iteration = %d\tmax_devia=%f\n',iteration,max_Devia);
% P2fine = P2;
% end
% len = length(P2);
% dp = 0.15;
% neighbor = 30;
% affmat = fitgeotrans(P1,P2,'projective');
% P3 = [P1 ones(len,1)] * affmat.T;
% P3 = P3(:,1:2) ./ P3(:,3);
% P2fine = P2;
% Origin_P2error = sqrt(sum((P2' - P3').^2));
% [~,ind] = sort(Origin_P2error);
% ref_match = ind(1:4);
% ind(1:4)=[];
% maxerror = 100;
% for i = 1:len-4
%     P2temp = P2fine([ref_match ind(i)],:); % P2fine is updated
%     for u = -neighbor:neighbor
%         du = dp * u;
%         for v = -neighbor:neighbor
%             dv = dp * v;
%             P2temp(end,1:2) = P2fine(ind(i),1:2) + [du dv]; % add subpixel deviation
%             affmat = fitgeotrans(P1([ref_match ind(i)],:),P2temp,'projective');
%             P3fine = [P1([ref_match ind(i)],:) ones(5,1)] * affmat.T;
%             P3fine = P3fine(:,1:2) ./ P3fine(:,3); % compute updated error of i(th) point
%             Delta_P2 = sum(sum((P3fine - P2temp)'.^2));
%             if Delta_P2 <= maxerror
%                 maxerror = Delta_P2;
%                 P2fine(ind(i),:) = P2temp(end,:); % update i(th) element of P2fine
%             end
%         end
%     end
% end
