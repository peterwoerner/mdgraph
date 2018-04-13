for i = 1:15
    r = Lattice2D(1, [i*10+20, i*10+20], 'sq2');
    r = GenerateHole(r, [i*10+20/2-0.25,i*10+20/2-0.25], 5);
    [natom(i), ~] = size(r)
    q = CreateChargeList(r);
    tic
    [Af, time1i, time2i] = combinedCoulombAlgo(r, 5, 1, q);
    time(i) = toc
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

plot(natom, natom*mean(time./(natom.*natom)),'--b', 'LineWidth', 3)
loglog(natom,time1, '-o', 'LineWidth', 3)
loglog(natom,time2, '-o', 'LineWidth', 3)
legend('Time Experiment', 'Quadratic Fit', 'Partition', 'Force Calc')

function [Af, time1, time2] = combinedCoulombAlgo(r, rc, k, q)
xmin = min(r(:,1));
ymin = min(r(:,2));

[rlength, ~] = size(r);
tic
kx = ceil((r(:,1)+xmin)/rc+eps);
ky = ceil((r(:,2)+ymin)/rc+eps);
%[min(min(r(:,1))), min(min(r(:,2))]
rpart{max(max(kx)), max(max(ky))} = [];
indices{max(max(kx)), max(max(ky))} = [];
%% Step 1 Bin Atoms
for j = 1:rlength
    rpart{kx(j),ky(j)} = [rpart{kx(j),ky(j)}; r(j,:)];
    indices{kx(j),ky(j)} = [indices{kx(j),ky(j)}; j];
end
time1 = toc
Ar = sparse(rlength,rlength);
Af = sparse(rlength,rlength);
nx = max(max(kx));
ny = max(max(ky));
n = nx*ny;

%For each bin and adjacent bin
tic
for ix = 1:nx
  for iy = 1:ny  %loop over all bins
    kx = ix;
    ky = iy;
    %indices{kx,ky}
    plength = length(indices{kx,ky});
    if ix == nx
        k1max = 0;
    else
        k1max = 1;
    end
    if iy == ny
        k2max = 0;
    else
        k2max = 1;
    end
        
        
    for k1 = 0:k1max     %loop over adjacent bins    
      for k2 = 0:k2max
        plength2 = length(indices{kx+k1,ky+k2});
        Arsub = zeros(plength2, plength);
        Afsub = zeros(plength2, plength);
        for j = 1:plength
            index=indices{kx,ky}(j);      
            Arsub(:,j) = sqrt((r(index,1) - r(indices{kx+k1,ky+k2}(:),1)).^2 + (r(index,2) - r(indices{kx+k1,ky+k2}(:),2)).^2);           
            Afsub(:,j) = k*q(index)*q(indices{kx+k1,ky+k2}(:))./(Arsub(:,j).^2+eps);
        end
        %Afsub = abs(Afsub + transpose(Afsub));
        %Arsub = Arsub + transpose(Arsub);
        Aff = zeros(max(plength,plength2));
        [xx,yy] = size(Afsub);
        Aff(1:xx, 1:yy) = Afsub;
        %Aff = sparsify_spectral(abs(Aff),1); %sparsify
        Afsub = Aff(1:xx, 1:yy);
        Af(indices{kx+k1,ky+k2},indices{kx,ky}) = Afsub; %transfer to true matrix
        Ar(indices{kx+k1,ky+k2},indices{kx,ky}) = Arsub; %transfer to true matrix
      end
    end
  end    
end
time2 = toc
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

