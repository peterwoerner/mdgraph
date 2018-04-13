for i = 1:12
    r = Lattice2D(1, [i*5+20, i*5+20], 'sq2');
    r = GenerateHole(r, [i*5+20/2-0.25,i*5+20/2-0.25], 5);
    [natom(i), ~] = size(r)
    q = abs(CreateChargeList(r));
    tic
    [~,Af, ~] = getACoulomb(r, q, 1);
    clear r q
    time(i) = toc
    [Af, time1i, time2i] = sparsify_spectral(Af,1);
    time1(i) = time1i
    time2(i) = time2i
end
loglog(natom,time, '-o', 'LineWidth', 3)
xlabel('Number of Atoms')
ylabel('Computational Time (s)')
hold on
% plot(natom, natom*mean(time./natom),'--r', 'LineWidth', 3)
% legend('Time Experiment', 'Linear Fit')
% %vars =  transpose(time)\[ones(15,1) transpose(natom.*log(natom))]
% plot(natom, natom.*log(natom)*mean(time./(natom.*log(natom))),'--k', 'LineWidth', 3)
% legend('Time Experiment', 'Linear Fit', 'Linear Log Fit')

plot(natom, natom.^3*mean(time1./(natom.^3)),'--b', 'LineWidth', 3)
loglog(natom,time1, '-o', 'LineWidth', 3)
loglog(natom,time2, '-o', 'LineWidth', 3)
legend('Time Experiment', 'Linear Fit', 'Partition', 'Force Calc')


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

