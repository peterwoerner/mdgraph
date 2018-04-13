function [] = MD2D( varargin)%[r, v]
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%Units discussion: Force: ev/Angstrom, position is in angstroms (10^-10 m),
%
coulombMethod = 'sparsify'
saveData = false;
saveDir = '/backupDisk1/mdtestHoleLong3';
mkdir(saveDir)
mkdir(strcat(saveDir,'/control'))
mkdir(strcat(saveDir,'/sparsify'))
mkdir(strcat(saveDir,'/threshold'))
mkdir(strcat(saveDir,'/combine'))
mkdir(strcat(saveDir,'/snapshot'))
mkdir(strcat(saveDir,'/final'))
saveFrequency = 1;
%updateFrequencies  = [0.5];%[-1,0.5,1,0,2.0,5];%[-1,0.5,1,0,2.0,5]%,10,30,100];
set(0, 'DefaultAxesFontSize', 20)
sparsificationFactor = 1;
coulombCut = 15.0;
graph = false;
tf =2000.0e-2; %final time
m = 1; %mass z%in g/mol
epsilon = 1;
sigma = 1;
omega0 = sqrt(72*epsilon/m/sigma^2/2^(1/3)) %approximate frequency of vibrations
%tf = 1/omega0
h = 1e-4;
while h > 1/omega0/10
    h = h/10;
end
t = 0:h:tf;
[ ~, Nt] = size(t)
Nt = Nt;
clear t
%mex forceCalcC.c
%mex forceCalcLJ.c
for jjjjj = 1:4
    if jjjjj == 1
        coulombMethod = 'sparsify'
        %coulombMethod = 'threshold'
    elseif jjjjj == 2
        coulombMethod = 'threshold'
        coulombCut = 15.0;
    elseif jjjjj == 3
        coulombMethod = 'combine'
        coulombCut = 15.0;
    else
        coulombMethod = 'threshold'
        coulombCut = 1e5;
    end
    for iiiii = 1:1
        strainRate = 1e-11/h;
        %% Initialize Run
        clear r
        global r v coulombList LJlist
        tic
        updateFrequency = pi;%ceil(0.5/omega0/h);
        %forceCF = 1.602*6.022e3; %Force conversion factor between eV/(g/mol) and \AA^2/ps^2 will give force/mass in \AA/ps^2, note that force will have ev/\AA
        %keConversion = 1/10^11/6.02/1.602; %conversion from g/mol*ps^2/\AA^2 to electron volt per particle
        forceCF = 1;
        keConversion = 1;

        xdim = 30*iiiii;
        ydim = 30*iiiii;
        r = Lattice2D(1, [xdim, ydim], 'sq2');
        aspectRatio = 1;
        b = 3.0/sqrt(aspectRatio);
        a = b*aspectRatio;
        r = GenerateEllipse(r, [0,ydim/2], a, b);
        topAtoms = find(r(:,2) > ydim - 1);
        bottomAtoms = find(r(:,2) < 1);
        sideAtoms = find(r(:,1) < 1);
        
        %r = [0, 0; 0.9, 0]
        r = r*1.56;
        rBottom = r(bottomAtoms,:);
        rTop = r(topAtoms,:);
        rSide = r(sideAtoms,:);
        % r = r + normrnd(0,.01, size(r))
        [Na, ~] = size(r)
        v = initvMag(r, 0);
        if graph == true
          figure
          scatter(r(:,1), r(:,2))
          figure
          scatter(r(:,1), r(:,2))
        end
        q = abs(CreateChargeList(r));
        %[rij, dir] = distancematrix(r,r);

        if saveData == true
           j = 0;
           tt(j+1) = 0;
           saveDataFunc()
        end

        [l,n] = size(r);

        if updateFrequency == 0
            LJlist = ones(Na)-eye(Na);
            coulombList = LJlist;
        else
            updatePairList();
        end



        paramLJ.epsilon = epsilon;
        paramLJ.sigma = sigma;
        paramLJ.type = 'lj';
        paramC.type = 'coulomb';
        paramC.k = 1e-2;
        paramC.q = q;

        [aljx, aljy] = forceCalcLJ(r, full(LJlist), epsilon, sigma);
        [acx, acy] = forceCalcC(r, full(coulombList), q, paramC.k);
        f = transpose([sum(aljx); sum(aljy)]) + transpose([sum(acx); sum(acy)]);

        %PV(1) = ke(1)/m/2+ 48*sum(sum(rij.*((rij+eps).^-13 + (rij+eps).^-7)))/2 + sum(sum(rij.*Af))/2;


        %% Loop time steps
         for i = 1:Nt+1
             %rall{i} = r;
             %vall{i} = v;
             a = forceCF*f./m;
             r = r + v*h + 1/2*a*h^2;
             if rem(i,updateFrequency) == 0
                updatePairList(); 
             end
             [aljx, aljy] = forceCalcLJ(r, full(LJlist), epsilon, sigma);
             [acx, acy] = forceCalcC(r, full(coulombList), q, paramC.k);
             f2 = transpose([sum(aljx); sum(aljy)]) + transpose([sum(acx); sum(acy)]);

             if rem(i,updateFrequency) == 0
                updatePairList(); 
             end

             a2 = forceCF*f2./m;
             v = v + h*(a+a2)/2;
             if saveData == true && rem(i,saveFrequency) == 0
                tt(j+1) = i*h;
                saveDataFunc()
             end

             applyDirichlet()
             f = f2;

             if rem(i,250) == 0
               if jjjjj == 4
                   save(strcat(saveDir,'/control/Savecontrol',num2str(i)))
               else
                  
                   save(strcat(saveDir,'/snapshot/Save',coulombMethod,num2str(i)))
               end
             end
         end
    %  tt = t

%% Post processing
    if graph == true
    te= ke + U; 
    hold on
    scatter(r(:,1), r(:,2), 'rx')
    figure 
    scatter(r(:,1), r(:,2), 'rx')
    end
    if graph == true && saveData == true
       figure
       set(gca, 'fontsize', 20)
       hold on
       subplot(3,1,1)
       plot(tt*h*omega0*saveFrequency,ke, 'r-o')
       ylabel('KE')
       subplot(3,1,2)
       plot(tt*h*omega0*saveFrequency,U, 'g-o')
       ylabel('potentialEnergy')
       subplot(3,1,3)
       plot(tt*h*omega0*saveFrequency,te, 'b-o')
       ylabel('KE + PE')
       figure
       plot(tt*h*omega0*saveFrequency, PV)
    %   figure
    %   plot(trace1(:,1), trace1(:,2), 'x')
    %    figure
    %    plot(tt, rij, 'b-o')
    end


    %[el, ~] = potentialEnergy(r,epsilon);

    %if graph == true
    %    figure
    %    hist(el)
    %end
    [sizes(iiiii),~] = size(r);
    times{jjjjj}(iiiii) = toc;
    %edges{jjjjj}(iiiii) = nnz();
    toc
    save(datesave(strcat(saveDir, '/Save',coulombMethod,num2str(Na))))
end
end
sizes
times{1}/Nt
times{2}/Nt
times{3}/Nt
figure
for i = 1:3
    hold on
    plot(sizes,times{i}/Nt)
end
% for i = 1:3
%     hold on
%     plot(sizes,edges{i}/Nt)
% end
%save('./mdtest2/Benchmark1')
    function applyDirichlet()
        strain = strainRate*i;
        r(bottomAtoms,2) = rBottom(:,2);
        r(sideAtoms,1) = rSide(:,1);
        r(topAtoms,:) = [rTop(:,1) rTop(:,2)+rTop(:,2)*strain];
    end

    function  updatePairList()
      if strcmp(coulombMethod, 'sparsify')
         [Ar, Af, ~] = getACoulomb(r, q, 1); %1 = k is the coefficient for F = kq_1,q_2/r^2, U = kq_1q_2/r
         Af = abs(Af);
         Afp = Af;
         Afn = Af;
         for jj = 1:Na
            for k = jj+1:Na
                if sign(q(jj)) == sign(q(k))
                    Afp(jj,k) = 0;
                    Afp(k,jj) = 0;
                else
                    Afn(jj,k) = 0;
                    Afn(k,jj) = 0;
                end
            end
         end
         Fsp = sparsify_spectral(Afp, sparsificationFactor);
         Fsn = sparsify_spectral(Afn, sparsificationFactor);
         Sp = Fsp./(Afp+eps(1));
         Sn = Fsn./(Afn+eps(1));
         coulombList = Sp + Sn;

      elseif strcmp(coulombMethod, 'threshold')
         Ar = distanceMatrixMex(r,r);
         coulombList = +(Ar<coulombCut.*(Ar > 0));
      elseif strcmp(coulombMethod, 'combine')
        Ar = distanceMatrixMex(r,r);
        Ar = Ar.*(+(Ar<coulombCut));  %test speed with and without this line
        Af = q*transpose(q)./(Ar.*Ar+eps);
        %[~, Af, Au] = getACoulomb(r, q, 1); %k is the coefficient for F = kq_1,q_2/r^2, U = kq_1q_2/r
        Af = abs(Af);
        %Ac = (rij < coulombCut);
        %Af = Af.*Ac;
        Afp = Af;
        Afn = Af;
        [n, ~] = size(Af);
        for ii = 1:Na
            for jj = ii+1:Na
                if sign(q(ii)) == sign(q(jj))
                    Afp(ii,jj) = 0;
                    Afp(jj,ii) = 0;
                else
                    Afn(ii,jj) = 0;
                    Afn(jj,ii) = 0;
                end
            end
        end
        %Fsp = sparsify_spectral(Afp, sparsificationFactor);
        Fsn = sparsify_spectral(Afn, sparsificationFactor);

        Sp = 0*Fsn;%Fsp./(Afp+eps(1));
        Sn = Fsn./(Afn+eps(1));
        coulombList = Sp + Sn;
      end
      LJlist = (+(Ar<2.5).*(Ar > 0));
    end

    function saveDataFunc()
       ke(j+1) = sum(m/2*(v(:,1).*v(:,1) + v(:,2).*v(:,2))); %kinetic energy 
       U(j+1) = U(r, epsilon, sigma);
       [~, Af, Au] = getACoulomb(r, q, 1);
       U(j+1) = U(j+1) + sum(sum(Au)); %Energy in epsilon electron volts\
       LJContribution = 48*sum(sum(rij.*((rij+eps(1)).^-13 + (rij+eps(1)).^-7)))/2;
       coulombContribution =  sum(sum(rij.*Af))/2;
       PV(j+1) = ke(j+1)/m/2 + LJContribution + coulombContribution;
       j = j + 1;        
    end


end

function [v] = initv0(r)
v = 0*r;
end


function [v] = initvMag(r,mag)
[n, dim] = size(r);
angles = linspace(0, 2*pi, n+1);
angles = angles(1:n);

phi = angles(randperm(n));
theta = angles(randperm(n));

if dim == 2
    x = mag*cos(phi);
    y = mag*sin(phi);
    v = [transpose(x) transpose(y)];
elseif dim == 3
    x = mag*sin(theta)*cos(phi);
    y = mag*sin(theta)*sin(phi);
    z = mag*cos(theta);
    v = [transpose(x) transpose(y) transpose(z)];
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


function [count] = getEdgeCount()

global coulombList
[n,~] = size(coulombList);
count = 0;
for i = 1:n
    count = count + length(coulombList{i,1});
end


end
