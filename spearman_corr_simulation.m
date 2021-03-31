%% Spearman Correlation Simulation

h_mat = kron([1 0; 0 1], triu(ones(8)) .* ~eye(8));
h_mat(1:8,9:16) = -1;

TestMat = zeros(10000,5);
for i = 1:10000

    A = triu(rand(16)) .* ~eye(16);
    B = triu(rand(16)) .* ~eye(16);
    Dat_A = A(h_mat ~= 0);
    Dat_B = B(h_mat ~= 0);
    TR_A = tiedrank(Dat_A);
    TR_B = tiedrank(Dat_B);
    
    TestMat(i,1) = fast_corr(TR_A,TR_B);
    
    Dat_A_p1 = A(h_mat > 0);
    Dat_B_p1 = B(h_mat > 0);
    
    Dat_A_p2 = A(h_mat < 0);
    Dat_B_p2 = B(h_mat < 0);
    
    TestMat(i,2) = fast_corr(tiedrank(Dat_A_p1),tiedrank(Dat_B_p1));
    TestMat(i,3) = fast_corr(tiedrank(Dat_A_p2),tiedrank(Dat_B_p2));
    
    A_nan = A; A_nan(A == 0) = NaN;
    B_nan = B; B_nan(B == 0) = NaN;
    
    TR_A = tiedrank(A_nan);
    TR_B = tiedrank(B_nan);
    
    TestMat(i,4) = fast_corr(TR_A(h_mat > 0),TR_B(h_mat > 0));
    TestMat(i,5) = fast_corr(TR_A(h_mat < 0),TR_B(h_mat < 0));
    
end

figure
subplot(1,3,1)
hist(TestMat(:,1),25)
set(gca,'xlim',[-0.8 0.8])
subplot(1,3,2)
hist(TestMat(:,2:3),25)
set(gca,'xlim',[-0.8 0.8])
subplot(1,3,3)
hist(TestMat(:,4:5),25)
set(gca,'xlim',[-0.8 0.8])


