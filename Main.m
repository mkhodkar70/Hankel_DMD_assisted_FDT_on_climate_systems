clc
clear all
close all

nd = 15;           %% number of datasets in the ensemble
M = 4;             %% sampling freq 1:6h; 4:24h
P = 1;             %% \tau of DMD or HDMD (day)
method  = 3;       %% 1: DMD, 2: DMD, 3: Hankel-DMD
tau_inf = 300;     %% maximum lag
forcing_type = 1;  %% 1: tropical, 2: subtropical

if(method == 1)     %% optimal parameters for each method
    q = 1;
    r = 400;
    r_fdt = 200;      %% # of EOFs used for dimension reduction of FDT
elseif(method == 2)
    q = 1;
    r = 400;
elseif(method == 3)
    q = 5;           %% delay
    r = 500;         %% DMD truncation value 
end

%% Domain properties

N_p = 39;             %% # of pressure-wise grid points
N_y = 48;             %% # of latitudinal grid points
L_p = 10000;          %% pressure-wise domain length
L_y = 180;            %% latitudinal domain length   

%% Loading projection modes and setting the forcing function

load('Robust_Modes_3M/Modes_Lagday1_q1_r400.mat', 'y', 'pf');

m = nd*200000;        %% length of training set (day)
p = pf;
y = y(N_y+1:end);
clear pf

f_u = zeros(N_y, N_p);
f_T = zeros(N_y, N_p);
if(forcing_type == 1)
    for(i = 1:N_y)
        for(j = 1:N_p)
            f_T(i, j) = 0.2*exp(-(p(j)/100-300)^2/100^2 - y(i)^2/20^2);
        end
    end
elseif(forcing_type == 2)
    for(i = 1:N_y)
        for(j = 1:N_p)
            f_T(i, j) = 0.2*exp(-(p(j)/100-450)^2/125^2 - (y(i)-25)^2/15^2);
        end
    end
elseif(forcing_type == 3)
    for(i = 1:N_y)
        for(j = 1:N_p)
            f_T(i, j) = 0.2*exp(-(p(j)/100-500)^2/125^2 - (y(i)-55)^2/25^2);
        end
    end
elseif(forcing_type == 4)
    for(i = 1:N_y)
        for(j = 1:N_p)
            f_u(i, j) = -0.2*exp(-(p(j)/100-350)^2/125^2 - (y(i)-35)^2/25^2);
            f_T(i, j) = 0.25*exp(-(p(j)/100-650)^2/125^2 - (y(i)-55)^2/25^2);
        end
    end
end

f_u = reshape(f_u, N_y*N_p, 1);
f_T = reshape(f_T, N_y*N_p, 1);

%% Loading the ground truth

text = ['Test_cases/Response', num2str(forcing_type), 'p.mat'];
load(text, 'uXTave', 'TXTave');
R_true_u = 0.5*(uXTave(end/2+1:end, :) + uXTave(end/2:-1:1, :));
R_true_T = 0.5*(TXTave(end/2+1:end, :) + TXTave(end/2:-1:1, :));

text = ['Test_cases/Response', num2str(forcing_type), 'n.mat'];
load(text, 'uXTave', 'TXTave');
R_true_u = 0.5*(R_true_u - 0.5*(uXTave(end/2+1:end, :) + uXTave(end/2:-1:1, :)));
R_true_T = 0.5*(R_true_T - 0.5*(TXTave(end/2+1:end, :) + TXTave(end/2:-1:1, :)));

%% Calculating the modes for the first segment of the data

sdU = zeros(N_p, nd); 
sdT = zeros(N_p, nd);   %% read averaged standard deviation at different pressure levels 

if(method > 1)
    for ens=1:(nd/2+1)
        tic
        ens
        load(['Data/CTLrun' num2str(ens) '/ZonalWind.mat'], 'u4xDaily', 'y', 'pf')
        load(['Data/CTLrun' num2str(ens) '/Temp.mat'], 'T4xDaily')
        np=size(u4xDaily,1)/2*size(u4xDaily,2);
        ns=size(u4xDaily,3); %number of samples

        N = floor(ns/M);
        if(ens==1)
            XX=zeros(2*np*q,(N-q+1-P)*(ceil(nd/2)+1)*2,'single');
            YY=zeros(2*np*q,(N-q+1-P)*(ceil(nd/2)+1)*2,'single');
        end        

        uN=zeros(48,39,floor(ns/M),'single');
        uS=uN;
        TS=uN;
        TN=uN;
        for n=M:M:ns
            uN(:,:,n/M)= u4xDaily(49:end,:,n);
            TN(:,:,n/M)= T4xDaily(49:end,:,n);
            uS(:,:,n/M)= u4xDaily(48:-1:1,:,n);
            TS(:,:,n/M)= T4xDaily(48:-1:1,:,n);
        end
        clear u4xDaily T4xDaily
        days=size(uN,3);

        disp('Remove mean ...')
        uNave=mean(uN,3);
        uSave=mean(uS,3);   
        TNave=mean(TN,3);
        TSave=mean(TS,3);
        uNa=uN-uNave;
        uSa=uS-uSave;
        TNa=TN-TNave;
        TSa=TS-TSave;
        clear uN uS TN TS

        disp('std ...')
        sduN=zeros(48,39,'single');
        sduS=sduN;
        sdTS=sduN;
        sdTN=sduN;
        for j=1:48
            for k=1:39
                sduN(j,k)= std(squeeze(uNa(j,k,:)));
                sdTN(j,k)= std(squeeze(TNa(j,k,:)));
                sduS(j,k)= std(squeeze(uSa(j,k,:)));
                sdTS(j,k)= std(squeeze(TSa(j,k,:)));
            end
        end

        disp('Weighting ...')
        for k=1:39
            sdU(k,ens)=mean(0.5*squeeze(sduN(:,k)+sduS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
            sdT(k,ens)=mean(0.5*squeeze(sdTN(:,k)+sdTS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
        end
        for n=1:days 
            for j=1:48
                uNa(j,:,n)=squeeze(uNa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
                uSa(j,:,n)=squeeze(uSa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
                TNa(j,:,n)=squeeze(TNa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
                TSa(j,:,n)=squeeze(TSa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
            end
        end

        disp('Vectorizing ...')
        XNa=zeros(2*np,days);
        XSa=XNa;
        for n=1:days
            XNa(:,n)=[(reshape(uNa(:,:,n),np,1));(reshape(TNa(:,:,n),np,1))];
            XSa(:,n)=[(reshape(uSa(:,:,n),np,1));(reshape(TSa(:,:,n),np,1))];
        end
        clear uNa uSa TNa TSa

        HNa=zeros(2*np*q,N-q+1,'single');
        HSa=HNa;
        for j=1:N-q+1
            for i=1:q
                HNa((i-1)*2*np+1:i*2*np,j)=XNa(:,i+(j-1));
                HSa((i-1)*2*np+1:i*2*np,j)=XSa(:,i+(j-1));
            end
        end
        clear XNa XSa

        XN0=HNa(:,1:end-P);
        XNp=HNa(:,1+P:end);

        XS0=HSa(:,1:end-P);
        %size(XS0)
        XSp=HSa(:,1+P:end);
        clear HSa HNa

        XX(:,(ens-1)*2*(N-q+1-P)+1:(ens)*2*(N-q+1-P)) = [XN0 XS0];
        YY(:,(ens-1)*2*(N-q+1-P)+1:(ens)*2*(N-q+1-P)) = [XNp XSp];

        clear XN0 XS0 XNp XSp
        toc
    end

    disp('SVD ...')
    tic
    [U, S, V] = svd(XX, 'econ');
    disp('SVD done')
    clear XX
    U = U(:, 1:r);
    S = S(1:r, 1:r);
    V = V(:, 1:r);
    AHDMD = U'*YY*V/S;
    [EVec,EVal] = eig(AHDMD);
    disp('EIG done')
    G = YY*V/S;
    clear YY V
    for n=1:size(EVal,1)
        phi(:,n)=G*EVec(:,n)/EVal(n,n);
    end
    [phi, EVal] = dmd_sorted(phi, EVal);
    phi         = my_compress(phi, q);

    toc

    disp('Saving ...')
    save(['Modes_Lagday' num2str(P) '_q' num2str(q) '_r' num2str(r) '_1.mat'], ... 
        'AHDMD', 'r', 'q', 'M', 'y', 'pf', 'S', 'phi', 'sdU', 'sdT', '-v7.3')
    clear AHDMD S phi

%% Calculating the modes for the second segment of the data
% Note: there are two datasets overlap between the two segments

    for ens=1:ceil(nd/2)+1
        tic
        ens+ceil(nd/2)-1
        load(['Data/CTLrun' num2str(ens+ceil(nd/2)-1) '/ZonalWind.mat'], 'u4xDaily', 'y', 'pf')
        load(['Data/CTLrun' num2str(ens+ceil(nd/2)-1) '/Temp.mat'], 'T4xDaily')
        np=size(u4xDaily,1)/2*size(u4xDaily,2);
        ns=size(u4xDaily,3); %number of samples

        N = floor(ns/M);
        if(ens==1)
            XX=zeros(2*np*q,(N-q+1-P)*(ceil(nd/2)+1)*2,'single');
            YY=zeros(2*np*q,(N-q+1-P)*(ceil(nd/2)+1)*2,'single');
        end        

        uN=zeros(48,39,floor(ns/M),'single');
        uS=uN;
        TS=uN;
        TN=uN;
        for n=M:M:ns
            uN(:,:,n/M)= u4xDaily(49:end,:,n);
            TN(:,:,n/M)= T4xDaily(49:end,:,n);
            uS(:,:,n/M)= u4xDaily(48:-1:1,:,n);
            TS(:,:,n/M)= T4xDaily(48:-1:1,:,n);
        end
        clear u4xDaily T4xDaily
        days=size(uN,3);

        disp('Remove mean ...')
        uNave=mean(uN,3);
        uSave=mean(uS,3);   
        TNave=mean(TN,3);
        TSave=mean(TS,3);
        uNa=uN-uNave;
        uSa=uS-uSave;
        TNa=TN-TNave;
        TSa=TS-TSave;
        clear uN uS TN TS
        if(forcing_type == 2 && nd == 8 && q > 1)
            CT = 0.75;
        end
        disp('std ...')
        sduN=zeros(48,39,'single');
        sduS=sduN;
        sdTS=sduN;
        sdTN=sduN;
        for j=1:48
            for k=1:39
                sduN(j,k)= std(squeeze(uNa(j,k,:)));
                sdTN(j,k)= std(squeeze(TNa(j,k,:)));
                sduS(j,k)= std(squeeze(uSa(j,k,:)));
                sdTS(j,k)= std(squeeze(TSa(j,k,:)));
            end
        end

        disp('Weighting ...')
        for k=1:39
            sdU(k,ens+ceil(nd/2)-1)=mean(0.5*squeeze(sduN(:,k)+sduS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
            sdT(k,ens+ceil(nd/2)-1)=mean(0.5*squeeze(sdTN(:,k)+sdTS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
        end
        for n=1:days 
            for j=1:48
                uNa(j,:,n)=squeeze(uNa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
                uSa(j,:,n)=squeeze(uSa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
                TNa(j,:,n)=squeeze(TNa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
                TSa(j,:,n)=squeeze(TSa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
            end
        end

        disp('Vectorizing ...')
        XNa=zeros(2*np,days);
        XSa=XNa;
        for n=1:days
            XNa(:,n)=[(reshape(uNa(:,:,n),np,1));(reshape(TNa(:,:,n),np,1))];
            XSa(:,n)=[(reshape(uSa(:,:,n),np,1));(reshape(TSa(:,:,n),np,1))];
        end
        clear uNa uSa TNa TSa

        HNa=zeros(2*np*q,N-q+1,'single');
        HSa=HNa;
        for j=1:N-q+1
            for i=1:q
                HNa((i-1)*2*np+1:i*2*np,j)=XNa(:,i+(j-1));
                HSa((i-1)*2*np+1:i*2*np,j)=XSa(:,i+(j-1));
            end
        end
        clear XNa XSa

        XN0=HNa(:,1:end-P);
        XNp=HNa(:,1+P:end);

        XS0=HSa(:,1:end-P);
        %size(XS0)
        XSp=HSa(:,1+P:end);
        clear HSa HNa

        XX(:,(ens-1)*2*(N-q+1-P)+1:(ens)*2*(N-q+1-P)) = [XN0 XS0];
        YY(:,(ens-1)*2*(N-q+1-P)+1:(ens)*2*(N-q+1-P)) = [XNp XSp];

        clear XN0 XS0 XNp XSp
        toc
    end
    disp('SVD ...')
    tic
    [U, S, V] = svd(XX, 'econ');
    disp('SVD done')
    clear XX
    U = U(:, 1:r);
    S = S(1:r, 1:r);
    V = V(:, 1:r);
    AHDMD = U'*YY*V/S;
    [EVec,EVal] = eig(AHDMD);
    disp('EIG done')
    G = YY*V/S;
    clear YY V
    for n=1:size(EVal,1)
        phi(:,n)=G*EVec(:,n)/EVal(n,n);
    end
    [phi, EVal] = dmd_sorted(phi, EVal);
    phi         = my_compress(phi, q);

    toc

    disp('Saving ...')
    save(['Modes_Lagday' num2str(P) '_q' num2str(q) '_r' num2str(r) '_2.mat'], ... 
        'AHDMD', 'r', 'q', 'M', 'y', 'pf', 'S', 'phi', 'sdU', 'sdT', '-v7.3')
    clear AHDMD S phi
end

%% Choosing EOFs calculated from the whole data or calculating the robust modes if it is a (Hankel-)DMD-enhanced FDT

if(method == 1)
    for ens=1:nd
        tic
        ens
        load(['Data/CTLrun' num2str(ens) '/ZonalWind.mat'], 'u4xDaily', 'y', 'pf')
        load(['Data/CTLrun' num2str(ens) '/Temp.mat'], 'T4xDaily')
        np=size(u4xDaily,1)/2*size(u4xDaily,2);
        ns=size(u4xDaily,3); %number of samples

        N = floor(ns/M);
        if(ens==1)
            XX=zeros(2*np*q,(N-q+1-P)*nd*2,'single');
            YY=zeros(2*np*q,(N-q+1-P)*nd*2,'single');
        end        

        uN=zeros(48,39,floor(ns/M),'single');
        uS=uN;
        TS=uN;
        TN=uN;
        for n=M:M:ns
            uN(:,:,n/M)= u4xDaily(49:end,:,n);
            TN(:,:,n/M)= T4xDaily(49:end,:,n);
            uS(:,:,n/M)= u4xDaily(48:-1:1,:,n);
            TS(:,:,n/M)= T4xDaily(48:-1:1,:,n);
        end
        clear u4xDaily T4xDaily
        days=size(uN,3);

        disp('Remove mean ...')
        uNave=mean(uN,3);
        uSave=mean(uS,3);   
        TNave=mean(TN,3);
        TSave=mean(TS,3);
        uNa=uN-uNave;
        uSa=uS-uSave;
        TNa=TN-TNave;
        TSa=TS-TSave;
        clear uN uS TN TS

        disp('std ...')
        sduN=zeros(48,39,'single');
        sduS=sduN;
        sdTS=sduN;
        sdTN=sduN;
        for j=1:48
            for k=1:39
                sduN(j,k)= std(squeeze(uNa(j,k,:)));
                sdTN(j,k)= std(squeeze(TNa(j,k,:)));
                sduS(j,k)= std(squeeze(uSa(j,k,:)));
                sdTS(j,k)= std(squeeze(TSa(j,k,:)));
            end
        end
        disp('Weighting ...')
        for k=1:39
            sdU(k,ens)=mean(0.5*squeeze(sduN(:,k)+sduS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
            sdT(k,ens)=mean(0.5*squeeze(sdTN(:,k)+sdTS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
        end
        for n=1:days 
            for j=1:48
                uNa(j,:,n)=squeeze(uNa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
                uSa(j,:,n)=squeeze(uSa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
                TNa(j,:,n)=squeeze(TNa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
                TSa(j,:,n)=squeeze(TSa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
            end
        end

        disp('Vectorizing ...')
        XNa=zeros(2*np,days);
        XSa=XNa;
        for n=1:days
            XNa(:,n)=[(reshape(uNa(:,:,n),np,1));(reshape(TNa(:,:,n),np,1))];
            XSa(:,n)=[(reshape(uSa(:,:,n),np,1));(reshape(TSa(:,:,n),np,1))];
        end
        clear uNa uSa TNa TSa

        HNa=zeros(2*np*q,N-q+1,'single');
        HSa=HNa;
        for j=1:N-q+1
            for i=1:q
                HNa((i-1)*2*np+1:i*2*np,j)=XNa(:,i+(j-1));
                HSa((i-1)*2*np+1:i*2*np,j)=XSa(:,i+(j-1));
            end
        end
        clear XNa XSa

        XN0=HNa(:,1:end-P);
        XNp=HNa(:,1+P:end);

        XS0=HSa(:,1:end-P);
        %size(XS0)
        XSp=HSa(:,1+P:end);
        clear HSa HNa

        XX(:,(ens-1)*2*(N-q+1-P)+1:(ens)*2*(N-q+1-P)) = [XN0 XS0];
        YY(:,(ens-1)*2*(N-q+1-P)+1:(ens)*2*(N-q+1-P)) = [XNp XSp];

        clear XN0 XS0 XNp XSp
        toc
    end

    disp('SVD ...')
    tic
    [U, S, V] = svd(XX, 'econ');
    disp('SVD done')
    clear XX
    U = U(:, 1:r);
else
    text = ['Modes_Lagday', num2str(P), '_q', num2str(q), '_r', num2str(r), '_1.mat'];
    load(text, 'phi', 'y', 'pf');

    phi_u_1 = phi(1:N_y*N_p, :);
    phi_T_1 = phi(N_y*N_p+1:end, :);

    clear phi

    text = ['Modes_Lagday', num2str(P), '_q', num2str(q), '_r', num2str(r), '_2.mat'];
    load(text, 'phi');

    phi_u_2 = phi(1:N_y*N_p, :);
    phi_T_2 = phi(N_y*N_p+1:end, :);
    clear phi

    %% Robust Modes
    
    trsh = 0.90;   %% threshold for robustness of modes
    k    = 0;      %% # of robust modes
    for(i = 1:r)
        for(j = 1:r)
            pcc_1 = corrcoef(real(phi_u_1(:, i)), real(phi_u_2(:, j)));
            pcc_1 = pcc_1(1, 2);

            pcc_2 = corrcoef(real(phi_T_1(:, i)), real(phi_T_2(:, j)));
            pcc_2 = pcc_2(1, 2);

            if(abs(pcc_1) > trsh && abs(pcc_2) > trsh)
                k = k + 1; 
                index_1(k) = i;
                index_2(k) = j;
                phi_u_rob(:, k) = phi_u_1(:, i); 
                phi_T_rob(:, k) = phi_T_1(:, i);
                break
            end
        end
    end
end

%% Calculating the covariance matrices

clc

if(method == 1)
    B = U(:, 1:r_fdt);
else    
    B = [phi_u_rob; phi_T_rob];
end    
r_fdt = size(B, 2);

disp(['number of robust modes = ', num2str(r_fdt)])

C_0   = zeros(2*size(B,2), 2*size(B,2));
C_tau = zeros(2*size(B,2), 2*size(B,2), tau_inf);

for ens=1:nd
    tic
    disp(['Ensemble =' num2str(ens)])
    disp('Reading data ...')
    load(['Data/CTLrun' num2str(ens) '/ZonalWind.mat'], 'u4xDaily', 'y', 'pf')
    load(['Data/CTLrun' num2str(ens) '/Temp.mat'], 'T4xDaily')
    np=size(u4xDaily,1)/2*size(u4xDaily,2);
    ns=size(u4xDaily,3); %number of samples

    uN=u4xDaily(49:end,:,M:M:ns);
    uS=u4xDaily(48:-1:1,:,M:M:ns);
    TN=T4xDaily(49:end,:,M:M:ns);
    TS=T4xDaily(48:-1:1,:,M:M:ns);
    clear u4xDaily T4xDaily
    days=size(uN,3);
    
    disp('Remove mean ...')
    uNa=uN-mean(uN,3);
    clear uN
    uSa=uS-mean(uS,3);
    clear uS
    TNa=TN-mean(TN,3);
    clear TN
    TSa=TS-mean(TS,3);
    clear TS

    disp('Weighting ...')
    %for k=1:39
    %    sdU(k,ens)=mean(0.5*squeeze(sduN(:,k)+sduS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
    %    sdT(k,ens)=mean(0.5*squeeze(sdTN(:,k)+sdTS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
    %end
    for n=1:days
       for j=1:48
           uNa(j,:,n)=squeeze(uNa(j,:,n)).*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
           uSa(j,:,n)=squeeze(uSa(j,:,n)).*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
           TNa(j,:,n)=squeeze(TNa(j,:,n)).*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
           TSa(j,:,n)=squeeze(TSa(j,:,n)).*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
       end
    end

    disp('Xa ...')
    XNa=zeros(2*np,days);
    XSa=XNa;

    for n=1:days
        XNa(:,n)=[double(reshape(uNa(:,:,n),np,1));double(reshape(TNa(:,:,n),np,1))];
        XSa(:,n)=[double(reshape(uSa(:,:,n),np,1));double(reshape(TSa(:,:,n),np,1))];
    end
    clear uNa uSa TNa TSa

    disp('Regressing ...')
    xNa=zeros(2*size(B,2),days);
    xSa=xNa;
    for n=1:days
        if(mod(n,25000)==0)
          disp(n)
          toc
        end
	    xNa(:,n) = [squeeze(mvregress(B(1:end/2,:),squeeze(XNa(1:end/2,n))));squeeze(mvregress(B(end/2+1:end,:),squeeze(XNa(end/2+1:end,n))))];
        xSa(:,n) = [squeeze(mvregress(B(1:end/2,:),squeeze(XSa(1:end/2,n))));squeeze(mvregress(B(end/2+1:end,:),squeeze(XSa(end/2+1:end,n))))];
    end
    clear XNa XSa

    disp('Covariance matrices ...')
    C_0 = C_0 + 0.5*(xNa*xNa' + xSa*xSa');
    for m=1:tau_inf
        C_tau(:,:,m) = C_tau(:,:,m) + 0.5*((xNa(:,1+m:days))*(xNa(:,1:days-m))' + (xSa(:,1+m:days))*(xSa(:,1:days-m))');
    end
    clear xNa xSa
end
C_0   = (1/nd)*C_0;
C_tau = (1/nd)*C_tau;

%% Calculating LRFs (for all lags)

LRF = zeros(2*r_fdt, 2*r_fdt, tau_inf);
I   = eye(2*r_fdt);
 
for i=0:tau_inf
    if(i == 0 || i == tau_inf)
        coeff = 0.5;
    else
        coeff = 1.0;
    end
    if(i >= 1)
    	cov_integral = cov_integral + coeff*(M/4.0)*squeeze(C_tau(:, :, i));    
        L            = cov_integral/C_0;
        LRF(:, :, i) = -L\I;
    else
        cov_integral = coeff*(M/4.0)*C_0;
    end
end

disp('Linear Response Functions ...')

%% Calculating the response

clc

m = nd*200000;        %% length of training set (day)
p = pf;
clear pf

B_u = B(1:end/2,:);
B_T = B(end/2+1:end,:);

if(norm(f_T) > 0)
    f_T_r = mvregress(B_T, f_T);
else
    f_T_r = zeros(r_fdt, 1);
end
if(norm(f_u) > 0)
    f_u_r = mvregress(B_u, f_u);
else
    f_u_r = zeros(r_fdt, 1);
end

f_r   = [f_u_r; f_T_r];

Err_opt = 10;
alpha = 1:(M/4):tau_inf;
for(i = 1:length(alpha))

    M_fdt = squeeze(LRF(:, :, floor(alpha(i))));
    B     = double(B(:, 1:r_fdt));       

    R_fdt_r = -(M_fdt\f_r);     
    R_fdt_u_temp = real(B_u*R_fdt_r(1:end/2));        %% FDT zonal wind response
    R_fdt_T_temp = real(B_T*R_fdt_r(end/2+1:end));    %% FDT temperature response
    
    R_fdt_u_temp = reshape(R_fdt_u_temp, N_y, N_p);        %% FDT temperature response
    R_fdt_T_temp = reshape(R_fdt_T_temp, N_y, N_p);        %% FDT zonal wind response

    Err_u = norm(reshape(R_true_u, N_p*N_y, 1) - reshape(R_fdt_u_temp, N_p*N_y, 1))/norm(reshape(R_true_u, N_p*N_y, 1));
    Err_T = norm(reshape(R_true_T, N_p*N_y, 1) - reshape(R_fdt_T_temp, N_p*N_y, 1))/norm(reshape(R_true_T, N_p*N_y, 1));   %% Calculating errors
    
    if(Err_u < Err_opt)
        alpha_opt = alpha(i);
        Err_opt = Err_u;
    end
end

clear R_true_u_temp R_true_T_temp

f_T_rb = real(B_T*f_T_r);
f_T_rb = reshape(f_T_rb, N_y, N_p); 

M_fdt = squeeze(LRF(:, :, floor(alpha_opt)));

R_fdt_r = -(M_fdt\f_r);     
R_fdt_u = real(B_u*R_fdt_r(1:end/2));
R_fdt_T = real(B_T*R_fdt_r(end/2+1:end));

%% Plotting

f_u     = reshape(f_u, N_y, N_p); 
f_T     = reshape(f_T, N_y, N_p);            %% FDT Velocity response
R_fdt_u = reshape(R_fdt_u, N_y, N_p); 
R_fdt_T = reshape(R_fdt_T, N_y, N_p);        %% FDT Velocity response


y = y(N_y+1:end);
figure(1);
subplot(1, 2, 1);
contourf(y, p/100, f_u');
hold on
axis ij
colorbar();
axis('square');
caxis([-0.2 0.2]);
colormap(b2r(-0.2, 0.2));
set(gca, 'Fontsize', 12)
xlabel('$\mu$', 'interpreter', 'latex');
ylabel('$p \ (hPa)$', 'interpreter', 'latex');
title('Mechanical Forcing', 'Fontsize', 14);
subplot(1, 2, 2);
contourf(y, p/100, f_T');
hold on
axis ij
colorbar();
axis('square');
caxis([-0.2 0.2]);
colormap(b2r(-0.2, 0.2));
set(gca, 'Fontsize', 12)
xlabel('$\mu$', 'interpreter', 'latex');
ylabel('$p \ (hPa)$', 'interpreter', 'latex');
title('Thermal Forcing', 'Fontsize', 14);
savefig('force.mat')


Err_u = norm(reshape(R_true_u, N_p*N_y, 1) - reshape(R_fdt_u, N_p*N_y, 1))/norm(reshape(R_true_u, N_p*N_y, 1));
Err_T = norm(reshape(R_true_T, N_p*N_y, 1) - reshape(R_fdt_T, N_p*N_y, 1))/norm(reshape(R_true_T, N_p*N_y, 1));   %% Calculating errors

fprintf('Error in zonal-wind is %4.2f per cent \n', Err_u*100);
fprintf('Error in temperature is %4.2f per cent \n', Err_T*100);

figure(2);
subplot(2, 2, 1);
contourf(y, p/100, R_fdt_u', [-3:0.2:3]);
hold on
axis ij
colorbar();
caxis([-3 3]);
colormap(b2r(-3, 3));
axis('square');
set(gca, 'Fontsize', 12)
xlabel('$\mu$', 'interpreter', 'latex');
ylabel('$p \ (hPa)$', 'interpreter', 'latex');
title('Zonal-wind (FDT)', 'Fontsize', 14);

subplot(2, 2, 2);
contourf(y, p/100, R_true_u',[-3:0.2:3]);
hold on
axis ij
colorbar();
caxis([-3 3]);
colormap(b2r(-3, 3));
axis('square');
set(gca, 'Fontsize', 12)
xlabel('$\mu$', 'interpreter', 'latex');
ylabel('$p \ (hPa)$', 'interpreter', 'latex');
title('Zonal-wind (GCM)', 'Fontsize', 14);

subplot(2, 2, 3);
contourf(y, p/100, R_fdt_T', [-3:0.2:3]);
hold on
axis ij
colorbar();
caxis([-3 3]);
colormap(b2r(-3, 3));
axis('square');
set(gca, 'Fontsize', 12)
xlabel('$\mu$', 'interpreter', 'latex');
ylabel('$p \ (hPa)$', 'interpreter', 'latex');
title('Temperature (FDT)', 'Fontsize', 14);

subplot(2, 2, 4);
contourf(y, p/100, R_true_T', [-3:0.2:3]);
hold on
axis ij
colorbar();
caxis([-3 3]);
colormap(b2r(-3, 3));
axis('square');
set(gca, 'Fontsize', 12)
xlabel('$\mu$', 'interpreter', 'latex');
ylabel('$p \ (hPa)$', 'interpreter', 'latex');
title('Temperature (GCM)', 'Fontsize', 14);
savefig('Response.fig');

save('Responses.mat', 'R_fdt_T', 'R_true_T', 'R_fdt_u', 'R_true_u', 'y', 'p')
