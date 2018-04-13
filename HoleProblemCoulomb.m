function [ output_args ] = HolePlacements()
%
 set(0,'DefaultAxesFontSize',30)
 set(groot,'defaultLineLineWidth',2)
epsilon = 1;
sigma = 1;
Fmax = 48*(-1/((23/7)^(1/6))^13 + 1/2/((23/7)^(1/6))^7);
r = Lattice2D(1, [40, 40], 'sq2');
r = GenerateHole(r, [19.75,19.75], 5);
q2 = CreateChargeList2(r);
q = CreateChargeList(r);
q = abs(q);
r = 1.56*r; %Moves close to equilibrium position
% r = r + normrnd(0, .1, size(r));
% natoms = length(r(:,1))
% for i = 1:natoms
%     r = max(max(r))*rand(natoms,2);
% end
% figure
% scatter(r(:,1), r(:,2))
% axis equal
%r = [0.5*r(:,1), 2*r(:,2)];
%scatter(r(:,1), r(:,2))
[Ar,Af, Au] = getACoulomb(r, q, epsilon);
max(max(Af));
min(min(Af));
% input('Continue?')
%Af = Af.*signForceCoul(r);

%[force, ~] = forcecalcCoulomb(r, q, epsilon);
% size(force);
% fx = force(:,:,1);
% fy = force(:,:,2);
sum(sum(Au));
k1 = find(q>0);
k2 = find(q<0);
r;
q;
k1;
k2;
rp1 = r(k1,:);
rp2 = r(k2,:);
figure(96)
scatter(rp1(:,1),rp1(:,2), 100, 'b+')
axis equal
hold on
scatter(rp2(:,1),rp2(:,2), 100, 'r*')
legend('Positive charges', 'Negative Charges')
dir =  '/home/peter/Documents/mdgraph/2DHole/figures4/'
name = 'AtomPositions';
str = 'CrystalSeparate';
set(figure(96), 'Position', [204   272   763   585])
saveas(96, strcat(dir, name, str))
saveas(96, strcat(dir, name, str), 'epsc')
% j1 = find(q2>0);
% j2 = find(q2<0);
% rp1 = r(j1,:);
% rp2 = r(j2,:);
% figure
% scatter(rp1(:,1),rp1(:,2), 100, 'b+')
% axis equal
% hold on
% scatter(rp2(:,1),rp2(:,2), 100, 'r*')


%% Sparisfy Spectral Potential Energy
% sparsifyTest(-Au, Ar, r)
% 
% sparsifyTest(fx, Ar, r)
% 
% sparsifyTest(fy, Ar, r)
% 
% sparsifyTest(sqrt(fx.*fx + fy.*fy), Ar, r)

% min(min(Af))
% max(max(Af))
% min(min(sqrt(fx.*fx + fy.*fy)))
% max(max(sqrt(fx.*fx + fy.*fy)))
% % 
% norm(sqrt(fx.*fx + fy.*fy)-Af)
fmaxSeparateTest(Af, Ar, r, q)


end


function [results] = sparsifyTest(fx, Ar, r)
%Compares and plots 

Fmax = max(max(abs(fx)));
Fmin = min(min(fx));

epsilon(1) = 0;
figure
histogram(sum(fx));
hold on;
ylabel('Instances')
xlabel('Degree')
sgn = sign(fx);
fxa = abs(fx);
Px(1) = sum(sum(fx));
Cx(1) = sum(sum(fx));

min(sum(fx));
max(sum(fx));
[fx transpose(sum(fx))];
error(1) = 0;
edge(1) = nonzeroCount(fx);
edgec(1) = edge(1);
errorc(1) = 0;
for i = 1:10
    epsilon(i+1) = i/10;
    Fs = sparsify_spectral(fxa, epsilon(i+1));
    Fc = DistanceCut(fx, Ar, 10*(1-i/10)+ 0.5);
    Fsa = Fs.*sgn;
    Px(i+1) = sum(sum(Fsa));
    Cx(i+1) = sum(sum(Fc));
    edge(i+1) =  nonzeroCount(Fs);
    edgec(i+1) = nonzeroCount(Fc);
    edgee(i+1) = nnz(Fs);
    edgece(i+1) = nnz(Fc);
    Pxdeg = sum(Fsa);
    Cxdeg = sum(Fc);
    if mod(i,2) == 0
        histogram(Pxdeg)%,'DisplayStyle','stairs')
    end%,'DisplayStyle','stairs')
    
    error(i+1) = norm((Pxdeg - sum(fx))/Fmax)/length(Pxdeg);
    errorc(i+1) = norm((Cxdeg - sum(fx))/Fmax)/length(Cxdeg);
    length(Pxdeg);
    hold on
end
edgee = edge - edgee
edgece = edgec - edgece
%legend('\epsilon = 0', '\epsilon = 0.1','\epsilon = 0.2','\epsilon = 0.3','\epsilon = 0.4','\epsilon = 0.5','\epsilon = 0.6','\epsilon = 0.7','\epsilon = 0.8','\epsilon = 0.9','\epsilon = 1.0')
legend('\epsilon = 0','\epsilon = 0.2','\epsilon = 0.4','\epsilon = 0.6','\epsilon = 0.8','\epsilon = 1.0')
figure
plot(1 - edge/edge(1), Px)
hold on
plot(1 - edgec/edgec(1), Cx)
ylabel('Net Force')
xlabel('Percent Edges Removed')
legend('Spectral Sparsification', 'Cutting Sparsification')
figure
plot(1 - edge/edge(1), error)
hold on
plot(1 - edgec/edgec(1), errorc)
ylabel('Normalized Error')
xlabel('Percent Edges Removed')
legend('Spectral Sparsification', 'Cutting Sparsification')
Px;
% figure
% C2 = DistanceCut(fx, Ar, 2.5);
% weightedGraphPlot(C2, r, sum(C2));
% 
% fxs = sgn.*sparsify_spectral(fxa, 1.0);
% figure
% weightedGraphPlot(fx, r, sum(fx));
% figure
% weightedGraphPlot(fxs, r, sum(fxs));
% fxs5 = sgn.*sparsify_spectral(fxa, 0.5);
% figure
% weightedGraphPlot(fxs5, r, sum(fxs5));
end


function [n] = nonzeroCount(A)
%%Counts the number of nonzero entries in symmetric matrix A excluding the diagonal.
n = 0;
[nx, ny] = size(A);
for i = 1:nx
    for j = i+1:ny
        if A(i,j) ~= 0
            n = n+1;
        end
    end
end

end

function [q] = CreateChargeList(r)
[nx, ny] = size(r);
q = zeros(nx,1);
for i=1:nx
    if mod(floor(2*r(i,1)),2) == 0
        q(i) = 1;
    else
        q(i) = -1;
    end
end
end

function [q] = CreateChargeList2(r)
[nx, ny] = size(r);
q = zeros(nx,1);
for i=1:nx
    if mod(i,2) == 0
        q(i) = 1;
    else
        q(i) = -1;
    end
end
end

function [] = fmaxTest(Af, Ar, r, q)
directory = '/home/woates/Documents/mdgraph/2DHole/figures4/';
str = 'Random';
[~,dir] = distancematrix(r,r);
Fmax = max(max(abs(Af)));
sgn = signForceCoul(q);
fx = Af.*dir(:,:,1).*sgn;
fy = Af.*dir(:,:,2).*sgn;
epsilon(1) = 0;
% figure(3)
% histogram(sum(Af));
% figure(4)
% plot(sum(fx))
% figure(5)
% plot(sum(fx))
% hold on;
ylabel('Instances')
xlabel('Degree')
sgn = sign(Af);
fxa = abs(Af);
Px(1) = sum(sum(Af));
Cx(1) = sum(sum(Af));
fxserror(1) = 0;
fyserror(1) = 0;
fxcerror(1) = 0;
fycerror(1) = 0;
Fxss(1) = sum(sum(fx));
Fyss(1) = sum(sum(fy));
Fxcs(1) = sum(sum(fx));
Fycsf(1) = sum(sum(fy));
error(1) = 0;
edge(1) = nonzeroCount(Af);
edgec(1) = edge(1);
errorc(1) = 0;
distance = [40, 30, 25, 20, 15, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.5, 2, 1.5, 1.0];
for i = 1:20
    epsilon(i+1) = i/20;
    Fs = sparsify_spectral(Af, epsilon(i+1));
    Fc = DistanceCut(Af, Ar, distance(i));
    data(i).Fs = Fs;
    data(i).Fc = Fc;
    fxs = Fs.*dir(:,:,1).*sgn;
    fys = Fs.*dir(:,:,2).*sgn;
    fxc = Fc.*dir(:,:,1).*sgn;
    fyc = Fc.*dir(:,:,2).*sgn;

    fxserror(i+1) = sqrt(sum((sum(fxs) - sum(fx)).^2)/Fmax)/length(sum(fx));
    fyserror(i+1) = sqrt(sum((sum(fys) - sum(fy)).^2)/Fmax)/length(sum(fx));
    fxcerror(i+1) = sqrt(sum((sum(fxc) - sum(fx)).^2)/Fmax)/length(sum(fx));
    fycerror(i+1) = sqrt(sum((sum(fyc) - sum(fy)).^2)/Fmax)/length(sum(fx));
    
    fxsperror(i+1) = sqrt(sum((sum(abs(fxs)) - sum(abs(fx))).^2)/Fmax)/length(sum(fx));
    fysperror(i+1) = sqrt(sum((sum(abs(fys)) - sum(abs(fy))).^2)/Fmax)/length(sum(fx));
    fxcperror(i+1) = sqrt(sum((sum(abs(fxc)) - sum(abs(fx))).^2)/Fmax)/length(sum(fx));
    fycperror(i+1) = sqrt(sum((sum(abs(fyc)) - sum(abs(fy))).^2)/Fmax)/length(sum(fx));    
    
    Fxss(i+1) = sum(sum(fxs));
    Fyss(i+1) = sum(sum(fys));
    Fxcs(i+1) = sum(sum(fxc));
    Fycs(i+1) = sum(sum(fyc));
    
    Fsa = Fs.*sgn;
    Px(i+1) = sum(sum(Fsa));
    Cx(i+1) = sum(sum(Fc));
    edge(i+1) =  nonzeroCount(Fs);
    edgec(i+1) = nonzeroCount(Fc);
    edgee(i+1) = nnz(Fs)/2;
    edgece(i+1) = nnz(Fc)/2;
    Pxdeg = sum(Fsa);
    Cxdeg = sum(Fc);
%     if mod(i,4) == 0
%         figure(3)
%         hold on
%         histogram(Pxdeg)%,'DisplayStyle','stairs')
%         figure(4)
%         hold on
%         plot(sum(fxs))
%         figure(5)
%         hold on
%         plot(sum(fxc))
%     end%,'DisplayStyle','stairs')
    
    error(i+1) = norm((Pxdeg - sum(Af))/Fmax)/length(Pxdeg);
    errorc(i+1) = norm((Cxdeg - sum(Af))/Fmax)/length(Cxdeg);
    length(Pxdeg);
    hold on
end
edgee
edgee = edgee-edge
edgec
edgece = edgece-edgec
%legend('\epsilon = 0', '\epsilon = 0.1','\epsilon = 0.2','\epsilon = 0.3','\epsilon = 0.4','\epsilon = 0.5','\epsilon = 0.6','\epsilon = 0.7','\epsilon = 0.8','\epsilon = 0.9','\epsilon = 1.0')
% legend('\epsilon = 0','\epsilon = 0.2','\epsilon = 0.4','\epsilon = 0.6','\epsilon = 0.8','\epsilon = 1.0')
% 
% legend('\epsilon = 0','\epsilon = 0.2','\epsilon = 0.4','\epsilon = 0.6','\epsilon = 0.8','\epsilon = 1.0')
% figure
% histogram(fxs, 21)
% hold on
% histogram(fys,21)
% histogram(DistanceCut(Af, Ar, 2.5).*dir(:,:,1).*sgn, 21)
% histogram(DistanceCut(Af, Ar, 2.5).*dir(:,:,2).*sgn, 21)
% histogram(fx, 21)
% histogram(fy,21)
% legend('fxs', 'fys', 'fxc', 'fyc', 'fx', 'fy')
% % legend('fxs', 'fys', 'fc', 'f')
% xlabel('Degree')
% ylabel('Instances')
% length(distance)
% 
% figure
% hold on
% histogram(fys,21)
% %histogram(fxc, 10)
% histogram(fyc,21)
% %histogram(fx, 10)
% histogram(fy,21)
% %legend('fxs', 'fys', 'fxc', 'fyc', 'fx', 'fy')
% legend('fs', 'fc', 'f')
% xlabel('Degree')
% ylabel('Instances')
find 
figure(97)
plot(1 - edge/edge(1), Px)
hold on
plot(1 - edgec/edgec(1), Cx)
ylabel('Net Force')
xlabel('Fraction of Edges Removed')
legend('Spectral Sparsification', 'Cutting Sparsification')
name = 'FSum';
set(figure(97), 'Position', [204   272   763   585])
saveas(97, strcat(directory, name, str))
saveas(97, strcat(directory, name, str), 'epsc')

figure(95)
plot(1 - edge/edge(1), error)
hold on
plot(1 - edgec/edgec(1), errorc)
ylabel('Normalized Error')
xlabel('Fraction of Edges Removed')
legend('Spectral Sparsification', 'Cutting Sparsification')
name = 'FError';
set(figure(95), 'Position', [204   272   763   585])
saveas(95, strcat(directory, name, str))
saveas(95, strcat(directory, name, str), 'epsc')
Px;
% figure
% plot(1-edge/edge(1),Fxcs)
% hold on
% plot(1-edge/edge(1),Fycs)
% plot(1-edge/edge(1),Fxss)
% plot(1-edge/edge(1),Fycs)
% ylabel('Total force directional force')
% xlabel('Fraction of Edges Removed')
% legend('f_x sparsification', 'f_y sparsification', 'f_x cutting', 'f_y cutting');
figure(98)
plot(1-edge/edge(1),fxserror)
hold on
plot(1-edge/edge(1),fyserror)
plot(1-edgec/edgec(1),fxcerror)
plot(1-edgec/edgec(1),fycerror)
ylabel('Error')
xlabel('Fraction of Edges Removed')
legend('f_x sparsification', 'f_y sparsification', 'f_x cutting', 'f_y cutting');
name = 'Components';
set(figure(98), 'Position', [204   272   763   585])
saveas(98, strcat(directory, name, str))
saveas(98, strcat(directory, name, str), 'epsc')


figure(101)
name = 'ComponentsMagnitude';
plot(1-edge/edge(1),fxsperror)
hold on
plot(1-edge/edge(1),fysperror)
plot(1-edgec/edgec(1),fxcperror)
plot(1-edgec/edgec(1),fycperror)
ylabel('Error')
xlabel('Fraction of Edges Removed')
legend('f_x sparsification', 'f_y sparsification', 'f_x cutting', 'f_y cutting');
set(figure(101), 'Position', [204   272   763   585])
saveas(101, strcat(directory, name, str))
saveas(101, strcat(directory, name, str), 'epsc')

figure(99)
name = 'ComponentsAvg';
plot(1- edge/edge(1),sqrt(fxserror.^2+fyserror.^2))
hold on
plot(1-edgec/edgec(1),sqrt(fxcerror.^2 + fycerror.^2))
ylabel('Normalized Error')
xlabel('Fraction of Edges Removed')
legend('Spectral Sparsification', 'Cut Sparsification');
set(figure(99), 'Position', [204   272   763   585])
saveas(99, strcat(directory, name, str))
saveas(99, strcat(directory, name, str), 'epsc')

figure(102)
name = 'ComponentsMagAvg';
plot(1- edge/edge(1),sqrt(fxsperror.^2+fysperror.^2))
hold on
plot(1-edgec/edgec(1),sqrt(fxcperror.^2 + fycperror.^2))
ylabel('Normalized Error')
xlabel('Fraction of Edges Removed')
legend('Spectral Sparsification', 'Cut Sparsification');
set(figure(102), 'Position', [204   272   763   585])
saveas(102, strcat(directory, name, str))
saveas(102, strcat(directory, name, str), 'epsc')


figure(100)
C2 = DistanceCut(Af, Ar, 2.5);
weightedGraphPlot(C2, r, sum(C2));
axis equal
saveas(100, 'graphCoulombCut')
saveas(100, 'graphCoulombCut', 'epsc')
xlim([29 32])
ylim([26 30])
saveas(100, 'graphCoulombCutZoom', 'epsc')
%fs = sgn.*sparsify_spectral(Af, 1.0);
figure(200)
weightedGraphPlot(Af, r, sum(Af));
axis equal
saveas(200, 'graphCoulomb')
saveas(200, 'graphCoulomb', 'epsc')
xlim([29 32])
ylim([26 30])
saveas(200, 'graphCoulombZoom', 'epsc')
figure(300)
weightedGraphPlot(Fs, r, sum(Fs));
axis equal
saveas(300, 'graphCoulombSparse')
saveas(300, 'graphCoulombSparse', 'epsc')
xlim([29 32])
ylim([26 30])
saveas(300, 'graphCoulombSparseZoom', 'epsc')
save('HoleProblemCoulombRaPositive')
end

function [sgn] =signForceCoul(q)
sgn = q*transpose(q);

end


function [] = fmaxSeparateTest(Af, Ar, r, q)
directory = '/home/peter/Documents/mdgraph/2DHole/figures4/';
str = 'Crystal';
[~,dir] = distancematrix(r,r);
Afp = Af;
Afn = Af;
[n, ~] = size(Af)
for i = 1:n
    for j = i+1:n
        if sign(q(i)) == sign(q(j))
            Afp(i,j) = 0;
            Afp(j,i) = 0;
        else
            Afn(i,j) = 0;
            Afn(j,i) = 0;
        end
    end
end

Fmax = max(max(abs(Af)));
sgn = signForceCoul(q);
fx = Afp.*dir(:,:,1)-Afn.*dir(:,:,1);
fy = Afp.*dir(:,:,2)-Afn.*dir(:,:,2);
epsilon(1) = 0;
figure(3)
histogram(sum(Af),25,'Normalization','probability');
% figure(4)
% plot(sum(fx))
% figure(5)
% plot(sum(fx))
% hold on;
ylabel('Instances')
xlabel('Degree')
sgn = sign(Af);
fxa = abs(Af);
Px(1) = sum(sum(Af));
Cx(1) = sum(sum(Af));
fxserror(1) = 0;
fyserror(1) = 0;
fxcerror(1) = 0;
fycerror(1) = 0;
Fxss(1) = sum(sum(fx));
Fyss(1) = sum(sum(fy));
Fxcs(1) = sum(sum(fx));
Fycsf(1) = sum(sum(fy));

FSinferror(1) = norm(Af - Af, inf)/norm(Af, inf);
FCinferror(1) = norm(Af - Af, inf)/norm(Af, inf);
FSfroerror(1) = norm(Af - Af, inf)/norm(Af, 'fro');
FCfroerror(1) = norm(Af - Af, inf)/norm(Af, 'fro');
FS1error(1) = norm(Af - Af, inf)/norm(Af, 1);
FC1error(1) = norm(Af - Af, inf)/norm(Af, 1);
FS2error(1) = norm(Af - Af, inf)/norm(Af, 2);
FC2error(1) = norm(Af - Af, inf)/norm(Af, 2);

fxsinferror(1) = norm(fx - fx, inf)/norm(fx, inf);
fxcinferror(1) = norm(fx - fx, inf)/norm(fx, inf);
fxsfroerror(1) = norm(fx - fx, inf)/norm(fx, 'fro');
fxcfroerror(1) = norm(fx - fx, inf)/norm(fx, 'fro');
fxs1error(1) = norm(fx - fx, inf)/norm(fx, 1);
fxc1error(1) = norm(fx - fx, inf)/norm(fx, 1);
fxs2error(1) = norm(fx - fx, inf)/norm(fx, 2);
fxc2error(1) = norm(fx - fx, inf)/norm(fx, 2);

fysinferror(1) = norm(fy - fy, inf)/norm(fy, inf);
fycinferror(1) = norm(fy - fy, inf)/norm(fy, inf);
fysfroerror(1) = norm(fy - fy, inf)/norm(fy, 'fro');
fycfroerror(1) = norm(fy - fy, inf)/norm(fy, 'fro');
fys1error(1) = norm(fy - fy, inf)/norm(fy, 1);
fyc1error(1) = norm(fy - fy, inf)/norm(fy, 1);
fys2error(1) = norm(fy - fy, inf)/norm(fy, 2);
fyc2error(1) = norm(fy - fy, inf)/norm(fy, 2);

error(1) = 0;
edge(1) = nonzeroCount(Af);
edgec(1) = edge(1);
errorc(1) = 0;
edgea(1) = nonzeroCount(Afp)
edger(1) = nonzeroCount(Afn)
distance = [40, 30, 25, 20, 15, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.5, 2, 1.5, 1.0];
for i = 1:1
    epsilon(i+1) = 1;
    %Fsp = sparsify_spectral(Afp, epsilon(i+1));
    Fsp = zeros(size(Afp));
    Fsn = sparsify_spectral(Afn, epsilon(i+1));
    Fs = Fsp - Fsn;
    Fpc = Afp.*(Ar<15);
    Fnc = Afn.*(Ar<15);
    Fc = Fpc - Fnc;
    data(i).Fs = Fs;
    data(i).Fc = Fc;
    fxs = Fsp.*dir(:,:,1)-Fsn.*dir(:,:,1);
    fys = Fsp.*dir(:,:,2)-Fsn.*dir(:,:,2);
    fxc = Fc.*dir(:,:,1).*sgn;
    fyc = Fc.*dir(:,:,2).*sgn;
    
    edgea(i+1) = nonzeroCount(Fpc);
    edger(i+1) = nonzeroCount(Fnc);
    errorPositive(i+1) = sum(sum(Fsp - Afp).^2)
    errorNegative(i+1) = sum(sum(Fsn - Afn).^2)
    [max(max(Afn)),max(max(Afp)), max(max(Fsp)), max(max(Fsn))]
    errorQS(i+1) = sum(sum(Fs - (Afp - Afn)).^2)
    errorQC(i+1) = sum(sum(Fc - (Afp - Afn)).^2)
    
    FSinferror(i+1) = norm(Fs - Af, inf)/norm(Af, inf);
    FCinferror(i+1) = norm(Fs - Af, inf)/norm(Af, inf);
    
    FSfroerror(i+1) = norm(Fs - Af, 'fro')/norm(Af, 'fro');
    FCfroerror(i+1) = norm(Fs - Af, 'fro')/norm(Af, 'fro');
    
    FS1error(i+1) = norm(Fs - Af, 1)/norm(Af, 1);
    FC1error(i+1) = norm(Fs - Af, 1)/norm(Af, 1);
    
    FS2error(i+1) = norm(Fs - Af, 2)/norm(Af, 2);
    FC2error(i+1) = norm(Fs - Af, 2)/norm(Af, 2);

    fxserror(i+1) = mean((sum(fxs) - sum(fx)).^2)
    fyserror(i+1) = mean((sum(fys) - sum(fy)).^2)
    fxcerror(i+1) = mean((sum(fxc) - sum(fx)).^2)
    fycerror(i+1) = mean((sum(fyc) - sum(fy)).^2)
    
    fxsperror(i+1) = sqrt(sum((sum(abs(fxs)) - sum(abs(fx))).^2)/Fmax)/length(sum(fx));
    fysperror(i+1) = sqrt(sum((sum(abs(fys)) - sum(abs(fy))).^2)/Fmax)/length(sum(fx));
    fxcperror(i+1) = sqrt(sum((sum(abs(fxc)) - sum(abs(fx))).^2)/Fmax)/length(sum(fx));
    fycperror(i+1) = sqrt(sum((sum(abs(fyc)) - sum(abs(fy))).^2)/Fmax)/length(sum(fx)); 
    
    fxsinferror(i+1) = norm(fxs - fx, inf)/norm(fx, inf);
    fxcinferror(i+1) = norm(fxc - fx, inf)/norm(fx, inf);
    
    fxsfroerror(i+1) = norm(fxs - fx, 'fro')/norm(fx, 'fro');
    fxcfroerror(i+1) = norm(fxs - fx, 'fro')/norm(fx, 'fro');
    
    fxs1error(i+1) = norm(fxs - fx, 1)/norm(fx, 1);
    fxc1error(i+1) = norm(fxs - fx, 1)/norm(fx, 1);
    
    fxs2error(i+1) = norm(fxs - fx, 2)/norm(fx, 2);
    fxc2error(i+1) = norm(fxs - fx, 2)/norm(fx, 2);
    
    fysinferror(i+1) = norm(fys - fy, inf)/norm(fy, inf);
    fycinferror(i+1) = norm(fyc - fy, inf)/norm(fy, inf);
    
    fysfroerror(i+1) = norm(fys - fy, 'fro')/norm(fy, 'fro');
    fycfroerror(i+1) = norm(fys - fy, 'fro')/norm(fy, 'fro');
    
    fys1error(i+1) = norm(fys - fy, 1)/norm(fy, 1);
    fyc1error(i+1) = norm(fys - fy, 1)/norm(fy, 1);
    
    fys2error(i+1) = norm(fys - fy, 2)/norm(fy, 2);
    fyc2error(i+1) = norm(fys - fy, 2)/norm(fy, 2);
    
    Fxss(i+1) = sum(sum(fxs));
    Fyss(i+1) = sum(sum(fys));
    Fxcs(i+1) = sum(sum(fxc));
    Fycs(i+1) = sum(sum(fyc));
    
    Fsa = Fs.*sgn;
    Px(i+1) = sum(sum(abs(Fs)));
    Cx(i+1) = sum(sum(abs(Fc)));
    edge(i+1) =  nonzeroCount(Fs);
    edgec(i+1) = nonzeroCount(Fc);
    edgee(i+1) = nnz(Fs);
    edgece(i+1) = nnz(Fc);
    Pxdeg = sum(abs(Fsa));
    Cxdeg = sum(Fc);
    if mod(i,5) == 0
        figure
        histogram(Pxdeg,25,'Normalization','probability')%,'DisplayStyle','stairs')
    end%,'DisplayStyle','stairs')
    
    error(i+1) = norm((Pxdeg - sum(Af))/Fmax)/length(Pxdeg);
    errorc(i+1) = norm((Cxdeg - sum(Af))/Fmax)/length(Cxdeg);
    length(Pxdeg);
    hold on
end
save('HoleProblemCoulombData')
edgee = edge-edgee
edgece = edgec-edgece
%legend('\epsilon = 0', '\epsilon = 0.1','\epsilon = 0.2','\epsilon = 0.3','\epsilon = 0.4','\epsilon = 0.5','\epsilon = 0.6','\epsilon = 0.7','\epsilon = 0.8','\epsilon = 0.9','\epsilon = 1.0')
% legend('\epsilon = 0','\epsilon = 0.2','\epsilon = 0.4','\epsilon = 0.6','\epsilon = 0.8','\epsilon = 1.0')
% 
% legend('\epsilon = 0','\epsilon = 0.2','\epsilon = 0.4','\epsilon = 0.6','\epsilon = 0.8','\epsilon = 1.0')
% figure
% histogram(fxs, 21)
% hold on
% histogram(fys,21)
% histogram(DistanceCut(Af, Ar, 2.5).*dir(:,:,1).*sgn, 21)
% histogram(DistanceCut(Af, Ar, 2.5).*dir(:,:,2).*sgn, 21)
% histogram(fx, 21)
% histogram(fy,21)
% legend('fxs', 'fys', 'fxc', 'fyc', 'fx', 'fy')
% % legend('fxs', 'fys', 'fc', 'f')
% xlabel('Degree')
% ylabel('Instances')
% length(distance)
% 
% figure
% hold on
% histogram(fys,21)
% %histogram(fxc, 10)
% histogram(fyc,21)
% %histogram(fx, 10)
% histogram(fy,21)
% %legend('fxs', 'fys', 'fxc', 'fyc', 'fx', 'fy')
% legend('fs', 'fc', 'f')
% xlabel('Degree')
% ylabel('Instances')

figure(88)
max(Px)
plot(1 - edge/edgec(1), errorQS/max(errorQC))
hold on
plot(1 - edgec/edge(1), errorQC/max(errorQC))
ylabel('Qincluded Force Error')
xlabel('Fraction of Edges Removed')
legend('Spectral Sparsification', 'Cutting Sparsification')
name = 'FSum';
set(figure(97), 'Position', [204   272   763   585])
saveas(97, strcat(directory, name, str))
saveas(97, strcat(directory, name, str), 'epsc')

figure(89)
plot(1 - edgea/edgea(1), errorPositive/sum(sum(Afp).^2))
hold on
plot(1 - edger/edger(1), errorNegative/sum(sum(Afn).^2))
%plot(1 - edger/edger(1), 2*errorNegative/sum(sum(Afn).^2))
ylabel('Net Force')
xlabel('Fraction of Edges Removed')
legend('Attractive Sparsification', 'Repulsive Sparsification', '2x Repulsive Sparsification')
name = 'FSum';
set(figure(97), 'Position', [204   272   763   585])
saveas(97, strcat(directory, name, str))
saveas(97, strcat(directory, name, str), 'epsc')

figure(97)
max(Px)
plot(1 - edge/edge(1), Px/max(Px))
hold on
plot(1 - edgec/edgec(1), Cx/max(Cx))
ylabel('Net Force')
xlabel('Fraction of Edges Removed')
legend('Spectral Sparsification', 'Cutting Sparsification')
name = 'FSum';
set(figure(97), 'Position', [204   272   763   585])
saveas(97, strcat(directory, name, str))
saveas(97, strcat(directory, name, str), 'epsc')

figure(95)
plot(1 - edge/edge(1), error)
hold on
plot(1 - edgec/edgec(1), errorc)
ylabel('Normalized Error')
xlabel('Fraction of Edges Removed')
legend('Spectral Sparsification', 'Cutting Sparsification')
name = 'FError';
set(figure(95), 'Position', [204   272   763   585])
saveas(95, strcat(directory, name, str))
saveas(95, strcat(directory, name, str), 'epsc')
Px;
% figure
% plot(1-edge/edge(1),Fxcs)
% hold on
% plot(1-edge/edge(1),Fycs)
% plot(1-edge/edge(1),Fxss)
% plot(1-edge/edge(1),Fycs)
% ylabel('Total force directional force')
% xlabel('Fraction of Edges Removed')
% legend('f_x sparsification', 'f_y sparsification', 'f_x cutting', 'f_y cutting');
figure(98)
plot(1-edge/edge(1),fxserror)
hold on
plot(1-edge/edge(1),fyserror)
plot(1-edgec/edgec(1),fxcerror)
plot(1-edgec/edgec(1),fycerror)
ylabel('Error')
xlabel('Fraction of Edges Removed')
legend('f_x sparsification', 'f_y sparsification', 'f_x cutting', 'f_y cutting');
name = 'Components';
set(figure(98), 'Position', [204   272   763   585])
saveas(98, strcat(directory, name, str))
saveas(98, strcat(directory, name, str), 'epsc')
% 
 
 figure(101)
 name = 'ComponentsMagnitude';
 plot(1-edge/edge(1),fxsperror)
 hold on
 plot(1-edge/edge(1),fysperror)
 plot(1-edgec/edgec(1),fxcperror)
 plot(1-edgec/edgec(1),fycperror)
 ylabel('Error')
 xlabel('Fraction of Edges Removed')
 legend('f_x sparsification', 'f_y sparsification', 'f_x cutting', 'f_y cutting');
 set(figure(101), 'Position', [204   272   763   585])
 saveas(101, strcat(directory, name, str))
 saveas(101, strcat(directory, name, str), 'epsc')
 
 figure(99)
 name = 'ComponentsAvg';
 plot(1- edge/edge(1),sqrt(fxserror.^2+fyserror.^2))
 hold on
 plot(1-edgec/edgec(1),sqrt(fxcerror.^2 + fycerror.^2))
 ylabel('Normalized Error')
 xlabel('Fraction of Edges Removed')
 legend('Spectral Sparsification', 'Cut Sparsification');
 set(figure(99), 'Position', [204   272   763   585])
 saveas(99, strcat(directory, name, str))
 saveas(99, strcat(directory, name, str), 'epsc')
 
 figure(102)
 name = 'ComponentsMagAvg';
 plot(1- edge/edge(1),sqrt(fxsperror.^2+fysperror.^2))
 hold on
% plot(1-edgec/edgec(1),sqrt(fxcperror.^2 + fycperror.^2))
% ylabel('Normalized Error')
% xlabel('Fraction of Edges Removed')
% legend('Spectral Sparsification', 'Cut Sparsification');
% set(figure(102), 'Position', [204   272   763   585])
% saveas(102, strcat(directory, name, str))
% saveas(102, strcat(directory, name, str), 'epsc')
% 
% 
% figure(102)
% name = 'ComponentsMagAvg';
% plot(1- edge/edge(1),errors)
% hold on
% plot(1-edgec/edgec(1),errorc)
% ylabel('Normalized Error')
% xlabel('Fraction of Edges Removed')
% legend('Spectral Sparsification', 'Cut Sparsification');
% set(figure(102), 'Position', [204   272   763   585])
% saveas(102, strcat(directory, name, str))
% saveas(102, strcat(directory, name, str), 'epsc')


% 
% figure(103)
% name = 'NetForceAdjacencyError';
% plot(1- edge/edge(1),FSinferror)
% hold on
% plot(1-edgec/edgec(1),FCinferror)
% plot(1- edge/edge(1),FSfroerror)
% plot(1- edgec/edgec(1),FCfroerror)
% plot(1- edge/edge(1),FS1error)
% plot(1- edgec/edgec(1),FC1error)
% plot(1- edge/edge(1),FS2error)
% plot(1- edgec/edgec(1),FC2error)
% ylabel('Relative Error')
% xlabel('Fraction of Edges Removed')
% legend('Spectral Sparsification (L_\infty)', 'Cut Sparsification (L_\infty)','Spectral Sparsification (L_f)', 'Cut Sparsification (L_f)','Spectral Sparsification (L_1)', 'Cut Sparsification (L_1)','Spectral Sparsification (L_2)', 'Cut Sparsification (L_2)');
% set(figure(103), 'Position', [204   272   763   585])
% saveas(103, strcat(directory, name, str))
% saveas(103, strcat(directory, name, str), 'epsc')
% 
% 
% figure(104)
% name = 'NetForceAdjacencyError';
% plot(1- edge/edge(1),fxsinferror,'Marker', 'hexagram')
% hold on
% plot(1-edgec/edgec(1),fxcinferror,'Marker', 'pentagram')
% plot(1- edge/edge(1),fxsfroerror,'Marker', 'x')
% plot(1- edgec/edgec(1),fxcfroerror,'Marker', 'o')
% plot(1- edge/edge(1),fxs1error,'Marker', '+')
% plot(1- edgec/edgec(1),fxc1error,'Marker', '*')
% plot(1- edge/edge(1),fxs2error,'Marker', 'diamond')
% plot(1- edgec/edgec(1),fxc2error,'Marker', 'square')
% ylabel('Relative Error')
% xlabel('Fraction of Edges Removed')
% legend('Spectral Sparsification (L_\infty)', 'Cut Sparsification (L_\infty)','Spectral Sparsification (L_f)', 'Cut Sparsification (L_f)','Spectral Sparsification (L_1)', 'Cut Sparsification (L_1)','Spectral Sparsification (L_2)', 'Cut Sparsification (L_2)');
% set(figure(104), 'Position', [204   272   763   585])
% saveas(104, strcat(directory, name, str))
% saveas(104, strcat(directory, name, str), 'epsc')

[rows, ~]  = find(abs(r(:,1)-30)<3);
[cols, ~]  = find(abs(r(:,2)-28)<3);
indices = intersect(rows,cols)

figure(100)
C2 = DistanceCut(Af, Ar, 2.5);
weightedGraphPlot(C2(indices,indices), r(indices,:), sum(C2));
axis equal
saveas(100, 'graphCoulombCut')
saveas(100, 'graphCoulombCut', 'epsc')
% xlim([29 32])
% ylim([26 30])
% saveas(100, 'graphCoulombCutZoom', 'epsc')
fs = sgn.*sparsify_spectral(Af, 1.0);
figure(200)
weightedGraphPlot(Af(indices,indices), r(indices,:), sum(Af));
axis equal
saveas(200, 'graphCoulomb')
saveas(200, 'graphCoulomb', 'epsc')
xlim([29 32])
ylim([26 30])
% saveas(200, 'graphCoulombZoom', 'epsc')
figure(300)
weightedGraphPlot(Fs(indices,indices), r(indices,:), sum(Fs));
axis equal
saveas(300, 'graphCoulombSparse')
saveas(300, 'graphCoulombSparse', 'epsc')
% xlim([29 32])
% ylim([26 30])
% saveas(300, 'graphCoulombSparseZoom', 'epsc')
% save('HoleProblemCoulombRaPositive')
end