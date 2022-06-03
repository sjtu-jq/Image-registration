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

