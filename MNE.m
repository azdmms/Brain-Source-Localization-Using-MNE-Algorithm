clear 
close all
clc 

load('C:\Users\Azadeh\OneDrive\Desktop\semester9\signal_processing_az_hajipour\Ex6\Lab 6\Lab 6\Codes and Data\Codes and Data\ElecPosXYZ.mat') ;

%Forward Matrix
ModelParams.R = [8 8.5 9.2] ; % Radius of diffetent layers
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3]; 
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;
%% problem 1
scatter3(LocMat(1 , :),LocMat(2 , :),LocMat(3 , :),15,'filled');
title('The location of dipoles (resolution = 1cm)') ;

%% problem 2
scatter3(LocMat(1 , :),LocMat(2 , :),LocMat(3 , :),15,'filled');
hold on;
Elecnames=cell(21,1);
for i = 1:21
    scatter3(9.2*ElecPos{i}.XYZ(1),9.2*ElecPos{i}.XYZ(2),9.2*ElecPos{i}.XYZ(3),20 , [0 0 0] ,'filled');
    text(9.2*ElecPos{i}.XYZ(1),9.2*ElecPos{i}.XYZ(2),9.2*ElecPos{i}.XYZ(3), ElecPos{i}.Name)
    Elecnames{i}=ElecPos{i}.Name;
end
title('Dipoles and electrodes location') ;
%% problem 3
scatter3(LocMat(1 , :),LocMat(2 , :),LocMat(3 , :),10,'filled');
hold on;
for i = 1:21
    scatter3(9.2*ElecPos{i}.XYZ(1),9.2*ElecPos{i}.XYZ(2),9.2*ElecPos{i}.XYZ(3),15 , [0 0 0] ,'filled');
    text(9.2*ElecPos{i}.XYZ(1),9.2*ElecPos{i}.XYZ(2),9.2*ElecPos{i}.XYZ(3), ElecPos{i}.Name)
end
loc_rand=LocMat(:,17);
x=loc_rand(1);
y=loc_rand(2);
z=loc_rand(3);
a=0.2;
u=x*a;
v=y*a;
w=z*a;
quiver3(x,y,z,u,v,w , 'LineWidth', 1.5 , 'ShowArrowHead','on');
title('Dipoles / Electrodes / Selected dipole')
%% problem4
load('C:\Users\Azadeh\OneDrive\Desktop\semester9\signal_processing_az_hajipour\Ex6\Lab 6\Lab 6\Codes and Data\Codes and Data\Interictal.mat') ;
Q=zeros(3951,10240);
Q(16*3+1,:)=x*Interictal(7,:);
Q(16*3+2,:)=y*Interictal(7,:);
Q(16*3+3,:)=z*Interictal(7,:);
M=GainMat*Q;
time = linspace(0 , length(Interictal(7,:))/256 , length(Interictal(7,:))) ;
figure ;
plot(time , Interictal(7,:))
xlabel('Time') ;
title('The spike activity of channel 7')
disp_eeg(M,35,256,Elecnames)
title('Electrodes Potential')
%% problem5
p1=[6.531 17.5 19.99 21.03 21.34 32.55 32.94 36.61];
p1=p1*256;
a=zeros(21,8,7);
for i=1:8
a(:,i,:) = M(: , p1(i)-3:p1(i)+3);
end
b = mean(mean(a , 2),3) ;
%
load('C:\Users\Azadeh\OneDrive\Desktop\semester9\signal_processing_az_hajipour\Ex6\Lab 6\Lab 6\Codes and Data\Codes and Data\ElecPatch.mat') ;
b=squeeze(b);
Display_Potential_3D(9.2,b);

%% problem6
G=GainMat;
alpha=0.5;
QMNE=G'*inv((G*G'+alpha*eye( 21 )))*M;
x1 = 10 ;
x2 = 30 ;
imagesc(time,linspace(x1,x2,x2-x1+1),QMNE(x1:x2,:));
title('Q_M_N_E (selected rows)')
colorbar ;
%% problem 7
apro = zeros(1317 , 10240) ;
for i = 1:1317
    apro(i , :) = sqrt(QMNE(3*(i-1)+1 , :).^2+QMNE(3*(i-1)+2 , :).^2+QMNE(3*(i-1)+3 , :).^2) ;
end
disp_eeg(apro(1:20 , :),5,256)
title('Approximated activity for each dipole')
%% problem 8
error = sqrt((LocMat(1 , 7)-LocMat(1 , 17))^2+(LocMat(2 , 7)-LocMat(2 , 17))^2+(LocMat(3 , 7)-LocMat(3 , 17))^2) ;
%% problem 9
scatter3(LocMat(1 , :),LocMat(2 , :),LocMat(3 , :),10,'filled');
hold on;
for i = 1:21
    scatter3(9.2*ElecPos{i}.XYZ(1),9.2*ElecPos{i}.XYZ(2),9.2*ElecPos{i}.XYZ(3),15 , [0 0 0] ,'filled');
    text(9.2*ElecPos{i}.XYZ(1),9.2*ElecPos{i}.XYZ(2),9.2*ElecPos{i}.XYZ(3), ElecPos{i}.Name)
end
di = 1000 ;
loc_rand=LocMat(:,di);
x=loc_rand(1);
y=loc_rand(2);
z=loc_rand(3);
a=0.2;
u=x*a;
v=y*a;
w=z*a;
quiver3(x,y,z,u,v,w , 'LineWidth', 1.5 , 'ShowArrowHead','on');
title('Dipoles / Electrodes / Selected dipole')
%%
Q=zeros(3951,10240);
Q((di-1)*3+1,:)=x*Interictal(7,:);
Q((di-1)*3+2,:)=y*Interictal(7,:);
Q((di-1)*3+3,:)=z*Interictal(7,:);
M=GainMat*Q;
time = linspace(0 , length(Interictal(7,:))/256 , length(Interictal(7,:))) ;
figure ;
plot(time , Interictal(7,:))
xlabel('Time') ;
title('The spike activity of channel 7')
disp_eeg(M,35,256,Elecnames)
title('Electrodes Potential')
%%
p1=[6.531 17.5 19.99 21.03 21.34 32.55 32.94 36.61];
p1=p1*256;
a=zeros(21,8,7);
for i=1:8
a(:,i,:) = M(: , p1(i)-3:p1(i)+3);
end
b = mean(mean(a , 2),3) ;
%
%load('C:\Users\Azadeh\OneDrive\Desktop\semester9\signal_processing_az_hajipour\Ex6\Lab 6\Lab 6\Codes and Data\Codes and Data\ElecPatch.mat') ;
b=squeeze(b);
Display_Potential_3D(9.2,b);
%%
G=GainMat;
alpha=0.5;
QMNE=G'*inv((G*G'+alpha*eye( 21 )))*M;
x1 = 900*3 ;
x2 = 1050*3 ;
imagesc(time,linspace(x1,x2,x2-x1+1),QMNE(x1:x2,:));
title('Q_M_N_E (selected rows)')
colorbar ;
%% 
apro = zeros(1317 , 10240) ;
for i = 1:1317
    apro(i , :) = sqrt(QMNE(3*(i-1)+1 , :).^2+QMNE(3*(i-1)+2 , :).^2+QMNE(3*(i-1)+3 , :).^2) ;
end
disp_eeg(apro(990:1010 , :),0.1,256)
title('Approximated activity for each dipole')
%% problem 8
error = sqrt((LocMat(1 , 7)-LocMat(1 , 17))^2+(LocMat(2 , 7)-LocMat(2 , 17))^2+(LocMat(3 , 7)-LocMat(3 , 17))^2) ;

%% problem 9
scatter3(LocMat(1 , :),LocMat(2 , :),LocMat(3 , :),10,'filled');
hold on;
for i = 1:21
    scatter3(9.2*ElecPos{i}.XYZ(1),9.2*ElecPos{i}.XYZ(2),9.2*ElecPos{i}.XYZ(3),15 , [0 0 0] ,'filled');
    text(9.2*ElecPos{i}.XYZ(1),9.2*ElecPos{i}.XYZ(2),9.2*ElecPos{i}.XYZ(3), ElecPos{i}.Name)
end
di = 800 ;
loc_rand=LocMat(:,di);
x=loc_rand(1);
y=loc_rand(2);
z=loc_rand(3);
a=0.2;
u=x*a;
v=y*a;
w=z*a;
quiver3(x,y,z,u,v,w , 'LineWidth', 1.5 , 'ShowArrowHead','on');
title('Dipoles / Electrodes / Selected dipole')
%%
Q=zeros(3951,10240);
Q((di-1)*3+1,:)=x*Interictal(7,:);
Q((di-1)*3+2,:)=y*Interictal(7,:);
Q((di-1)*3+3,:)=z*Interictal(7,:);
M=GainMat*Q;
time = linspace(0 , length(Interictal(7,:))/256 , length(Interictal(7,:))) ;
figure ;
plot(time , Interictal(7,:))
xlabel('Time') ;
title('The spike activity of channel 7')
disp_eeg(M,35,256,Elecnames)
title('Electrodes Potential')
%%
p1=[6.531 17.5 19.99 21.03 21.34 32.55 32.94 36.61];
p1=p1*256;
a=zeros(21,8,7);
for i=1:8
a(:,i,:) = M(: , p1(i)-3:p1(i)+3);
end
b = mean(mean(a , 2),3) ;
%
%load('C:\Users\Azadeh\OneDrive\Desktop\semester9\signal_processing_az_hajipour\Ex6\Lab 6\Lab 6\Codes and Data\Codes and Data\ElecPatch.mat') ;
b=squeeze(b);
Display_Potential_3D(9.2,b);
%%
G=GainMat;
alpha=0.5;
QMNE=G'*inv((G*G'+alpha*eye( 21 )))*M;
x1 = 790*3 ;
x2 = 810*3 ;
imagesc(time,linspace(x1,x2,x2-x1+1),QMNE(x1:x2,:));
title('Q_M_N_E (selected rows)')
colorbar ;
%% 
apro = zeros(1317 , 10240) ;
for i = 1:1317
    apro(i , :) = sqrt(QMNE(3*(i-1)+1 , :).^2+QMNE(3*(i-1)+2 , :).^2+QMNE(3*(i-1)+3 , :).^2) ;
end
disp_eeg(apro(790:810 , :),1,256)
title('Approximated activity for each dipole')
%% problem 8
error = sqrt((LocMat(1 , 7)-LocMat(1 , 17))^2+(LocMat(2 , 7)-LocMat(2 , 17))^2+(LocMat(3 , 7)-LocMat(3 , 17))^2) ;
%% functions
function [gridpoints,TransferMat] = ForwardModel_3shell(Resolution, ModelParams)

%Grid Locations
Radius = ModelParams.R(1) ;
[X,Y,Z] = meshgrid(-Radius:Resolution:Radius,-Radius:Resolution:Radius,-Radius*.1:Resolution:Radius) ;
X = reshape(X,1,[]) ;
Y = reshape(Y,1,[]) ;
Z = reshape(Z,1,[]) ;
grid = find(X.^2+Y.^2+Z.^2<Radius^2) ;
gridpoints = [X(grid); Y(grid); Z(grid)] ;
clear('X','Y','Z') ;
GridLocs_temp = reshape(gridpoints,1,3,[]) ;
GridLocs_temp = repmat(GridLocs_temp,[21,1,1]) ;

%Electrode Positions
load('C:\Users\Azadeh\OneDrive\Desktop\semester9\signal_processing_az_hajipour\Ex6\Lab 6\Lab 6\Codes and Data\Codes and Data\ElecPosXYZ.mat') ;
ElectrodePos = [] ;
for i=1:21
    ElectrodePos(i,:) = Radius*ElecPos{i}.XYZ ;
end
ElectrodePos = repmat(ElectrodePos,[1,1,size(GridLocs_temp,3)]) ;

TransferMat = zeros(size(GridLocs_temp)) ;

for shell=1:3
    GridLocs = ModelParams.Mu(shell)*GridLocs_temp ;
    
    C1 = 2*sum((ElectrodePos-GridLocs).*GridLocs,2)./(sum((ElectrodePos-GridLocs).^2,2)).^1.5 + ...
    1./sqrt(sum((ElectrodePos-GridLocs).^2,2))- ...
    1./sqrt(sum((ElectrodePos).^2,2)) ;
    C2 = 2./(sum((ElectrodePos-GridLocs).^2,2)).^1.5 + ...
    (sqrt(sum((ElectrodePos-GridLocs).^2,2))+ sqrt(sum((ElectrodePos).^2,2)))./...
    (sqrt(sum((ElectrodePos).^2,2)).*sqrt(sum((ElectrodePos-GridLocs).^2,2)).*...
    (sqrt(sum((ElectrodePos).^2,2)).*sqrt(sum((ElectrodePos-GridLocs).^2,2))+...
    sum((ElectrodePos).^2,2)-sum(GridLocs.*ElectrodePos,2))) ;
    h = (repmat(C1-C2.*sum(GridLocs.*ElectrodePos,2),[1,3,1]).*GridLocs + ...
        repmat(C2.*sum(GridLocs.^2,2),[1,3,1]).*ElectrodePos)./...
        (repmat(4*pi*ModelParams.Sigma(3)*sum(GridLocs.^2,2),[1,3,1])) ;
    TransferMat = TransferMat + ModelParams.Lambda(shell)*h ;
end

TransferMat = reshape(TransferMat,21,[]) ;

end


function t = disp_eeg(X,offset,feq,ElecName,titre)
% function t = disp_eeg(X,offset,feq,ElecName,titre)
%
% inputs
%     X: dynamics to display. (nbchannels x nbsamples) matrix
%     offset: offset between channels (default max(abs(X)))
%     feq: sapling frequency (default 1)
%     ElecName: cell array of electrode labels (default {S1,S2,...})
%     titre: title of the figure
%
% output
%     t: time vector
%
% G. Birot 2010-02


%% Check arguments
[N K] = size(X);

if nargin < 4
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
    titre = [];
end

if nargin < 5
    titre = [];
end

if isempty(feq)
    feq = 1;
end

if isempty(ElecName)
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
end

if isempty(offset)
    offset = max(abs(X(:)));
end


%% Build dynamic matrix with offset and time vector
X = X + repmat(offset*(0:-1:-(N-1))',1,K);
t = (1:K)/feq;
graduations = offset*(0:-1:-(N-1))';
shiftvec = N:-1:1;
Ysup = max(X(1,:)) + offset;
Yinf = min(X(end,:)) - offset;
% YLabels = cell(N+2) ElecName(shiftvec)

%% Display
figure1 = figure;
% a1 = axes('YAxisLocation','right');
a2 = axes('YTickLabel',ElecName(shiftvec),'YTick',graduations(shiftvec),'FontSize',7);
ylim([Yinf Ysup]);
box('on');
grid('on')
hold('all');
plot(t,X');
xlabel('Time (seconds)','FontSize',10);
ylabel('Channels','FontSize',10);
title(titre);
hold off
end

function Display_Potential_3D(Radius,ElecPot)

load('C:\Users\Azadeh\OneDrive\Desktop\semester9\signal_processing_az_hajipour\Ex6\Lab 6\Lab 6\Codes and Data\Codes and Data\ElecPosXYZ.mat') ;
load('C:\Users\Azadeh\OneDrive\Desktop\semester9\signal_processing_az_hajipour\Ex6\Lab 6\Lab 6\Codes and Data\Codes and Data\ElecPatch.mat') ;

%Nasion, Inion, Ears
Nz = Radius*[0.927+.25	 -0	-0.375] ;
Iz = Radius*[-0.906-.4	-1.11e-16 -0.423] ;

Ver = [] ;
for i=1:21
    Ver(i,:) = Radius*ElecPos{i}.XYZ ;
end
Face = ElecPatch ;

Cmat = reshape(ElecPot,[],1) ;
patch('Vertices',Ver,'Faces',Face,'FaceVertexCData',Cmat,'FaceColor','interp','FaceAlpha',1)
hold on
for i=1:21
    text(Ver(i,1) ,Ver(i,2), Ver(i,3), ElecPos{i}.Name) ;
end
text(Nz(1),Nz(2),Nz(3), 'Nasion','HorizontalAlignment','center','EdgeColor','r') ;
text(Iz(1),Iz(2),Iz(3), 'Inion','HorizontalAlignment','center','EdgeColor','r') ;
xlim([-1.5,1.5]*Radius) ;
ylim([-1.2,1.2]*Radius) ;
zlim([-.5,1.3]*Radius) ;
view(90,90) ;
colorbar
end