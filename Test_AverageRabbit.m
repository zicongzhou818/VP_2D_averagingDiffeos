clc
clear
close all
%% Reconstruct
N=101;h=1;

%% given phi
load Phi1_03rabbit
load Phi2_03rabbit
load cv_03rabbit
load jd_03rabbit

Phi1=Phi1_03rabbit;
Phi2=Phi2_03rabbit;

figure(1)
imshow(imresize(zeros(N)+255,[513 513]),'border', 'tight'); hold on
for i = 1:h:N
    plot(-Phi1(i,1:h:end),Phi2(i,1:h:end),'k-'); hold on
    plot(-Phi1(1:h:end,i),Phi2(1:h:end,i),'k-'); hold on;
end
axis([-N,-1,1,N]);
figure(2)
imshow(imresize(jd_03rabbit',[513 513]),'border', 'tight'); hold on
figure(3)
imshow(imresize(cv_03rabbit',[513 513]),'border', 'tight'); hold on

%% create 2 transformations
theta=pi/6;
cute=0.5;
X_new1=zeros(N);Y_new1=zeros(N);
X_new2=zeros(N);Y_new2=zeros(N);
h=1;
m=N;n=N;
x=1:n;
y=1:m;
[X,Y]=ndgrid(x,y);
for i=1:N 
    for j=1:N
      %% Rotate +
      [x_new, y_new]=cut_off_rotation(Phi1(i,j),Phi2(i,j),N,theta,cute);
      X_new1(i,j)=x_new;
      Y_new1(i,j)=y_new;
      %% Rotate -
      [x_new, y_new]=cut_off_rotation(Phi1(i,j),Phi2(i,j),N,-theta,cute);
      X_new2(i,j)=x_new;
      Y_new2(i,j)=y_new;
    end
end

figure(4)
imshow(imresize(zeros(N)+255,[513 513]),'border', 'tight'); hold on
for i = 1:h:N
    plot(-X_new1(i,1:h:end),Y_new1(i,1:h:end),'k-'); hold on
    plot(-X_new1(1:h:end,i),Y_new1(1:h:end,i),'k-'); hold on;
end
axis([-N,-1,1,N]);
figure(5)
imshow(imresize(zeros(N)+255,[513 513]),'border', 'tight'); hold on
for i = 1:h:N
    plot(-X_new2(i,1:h:end),Y_new2(i,1:h:end),'k-'); hold on
    plot(-X_new2(1:h:end,i),Y_new2(1:h:end,i),'k-'); hold on;
end
axis([-N,-1,1,N]);
[JD_new1,CV_new1]=compute_JD_and_Curl(X_new1,Y_new1,h);
[JD_new2,CV_new2]=compute_JD_and_Curl(X_new2,Y_new2,h);
figure(6)
imshow(imresize(JD_new1',[513 513]),'border', 'tight'); hold on
figure(7)
imshow(imresize(JD_new2',[513 513]),'border', 'tight'); hold on
figure(8)
imshow(imresize(CV_new1',[513 513]),'border', 'tight'); hold on
figure(9)
imshow(imresize(CV_new2',[513 513]),'border', 'tight'); hold on
%% reconstruct avg
JD_avg=nthroot(JD_new1.*JD_new2,2); 
CV_avg=(CV_new1+CV_new2)/2;
maaaxJD=max(max(JD_avg))
miiixJD=min(min(JD_avg))
maaaxCV=max(max(CV_avg))
miiixCV=min(min(CV_avg))

tic
[Avg1X,Avg1Y,~,~,~,~,~,~,Rot1]=PJDC_on_given_mesh2(JD_avg,CV_avg,N,X,Y);
toc

figure(10)
imshow(imresize(zeros(N)+255,[513 513]),'border', 'tight'); hold on
for i = 1:h:N
    plot(-Avg1X(i,1:h:end),Avg1Y(i,1:h:end),'r-'); hold on
    plot(-Avg1X(1:h:end,i),Avg1Y(1:h:end,i),'r-'); hold on;
end
axis([-N,-1,1,N]);

figure(11)
imshow(imresize(zeros(N)+255,[513 513]),'border', 'tight'); hold on
for i = 1:h:N
    plot(-Phi1(i,1:h:end),Phi2(i,1:h:end),'k-'); hold on
    plot(-Phi1(1:h:end,i),Phi2(1:h:end,i),'k-'); hold on;
end
for i = 1:h:N
    plot(-Avg1X(i,1:h:end),Avg1Y(i,1:h:end),'r-'); hold on
    plot(-Avg1X(1:h:end,i),Avg1Y(1:h:end,i),'r-'); hold on;
end
axis([-N,-1,1,N]);

[Avg1_JD,Avg1_CV]=compute_JD_and_Curl(Avg1X,Avg1Y,h);
maaaxJDdiff1=(max(max((Avg1_JD-jd_03rabbit))))
maaaxCVdiff1=(max(max((Avg1_CV-cv_03rabbit))))
maaaxdiff1=(max(max(((Phi1-Avg1X).^2+(Phi2-Avg1Y).^2).^0.5)))

%% average rabbit
tic
[Avg2X,Avg2Y,~,~,~,~,~,~,Rot2]=PJDC_on_given_mesh2fast(JD_avg,CV_avg,N,X,Y);
toc

figure(12)
imshow(imresize(zeros(N)+255,[513 513]),'border', 'tight'); hold on
for i = 1:h:N
    plot(-Avg2X(i,1:h:end),Avg2Y(i,1:h:end),'r-'); hold on
    plot(-Avg1X(1:h:end,i),Avg1Y(1:h:end,i),'r-'); hold on;
end
axis([-N,-1,1,N]);

figure(13)
imshow(imresize(zeros(N)+255,[513 513]),'border', 'tight'); hold on
for i = 1:h:N
    plot(-Phi1(i,1:h:end),Phi2(i,1:h:end),'k-'); hold on
    plot(-Phi1(1:h:end,i),Phi2(1:h:end,i),'k-'); hold on;
end
for i = 1:h:N
    plot(-Avg2X(i,1:h:end),Avg2Y(i,1:h:end),'r-'); hold on
    plot(-Avg2X(1:h:end,i),Avg2Y(1:h:end,i),'r-'); hold on;
end
axis([-N,-1,1,N]);

[Avg2_JD,Avg2_CV]=compute_JD_and_Curl(Avg2X,Avg2Y,h);
maaaxJDdiff2=(max(max((Avg2_JD-jd_03rabbit))))
maaaxCVdiff2=(max(max((Avg2_CV-cv_03rabbit))))
maaaxdiff2=(max(max(((Phi1-Avg2X).^2+(Phi2-Avg2Y).^2).^0.5)))

mm=length(Rot2);
figure(14)
plot(Rot1(1:mm,1),'b-'), ylabel('ratios'), title('ratio_{alg1} v.s. ratio_{alg2}'), hold on
plot(Rot2,'g-'), xlabel('Number of iterations'), legend('ratio_{alg1}','ratio_{alg2}'), grid on
a=axis;axis([-a(2)/50 a(2) -a(4)/50 a(4)]);

figure(15)
plot(log(Rot1(1:mm,1)),'b-'), ylabel('log of ratios'),title('ln(ratio_{alg1}) and ln(ratio_{alg2}) on iterations'), hold on
plot(log(Rot2),'g-'), xlabel('Number of iterations'), legend('ln(ratio_{alg1})','ln(ratio_{alg2})'), grid on
a=axis;axis([-a(2)/50 a(2) a(3) a(4)]);

figure(16)
line_fewer_markers(1:mm,Rot1(1:mm,1),20,'-bv','spacing','curve','linewidth',1,'markersize',6); hold on
line_fewer_markers(1:mm,Rot2,20,'-go','spacing','curve','linewidth',1,'markersize',6); grid on
legend('ratio_{alg1}','ratio_{alg2}')
xlabel('Number of iterations');ylabel('ratios');
a=axis;axis([-a(2)/50 a(2) -a(4)/50 a(4)]);


figure(17)
line_fewer_markers(1:mm,log(Rot1(1:mm,1)),20,'-bv','spacing','curve','linewidth',1,'markersize',6); hold on
line_fewer_markers(1:mm,log(Rot2),20,'-go','spacing','curve','linewidth',1,'markersize',6); grid on
legend('log ratio_{alg1}','log ratio_{alg2}')
xlabel('Number of iterations');ylabel('log ratios');
a=axis;axis([-a(2)/50 a(2) a(3) a(4)]);



% maaaxJD =
% 
%    1.808417515010706
% 
% 
% miiixJD =
% 
%    0.582625381728722
% 
% 
% maaaxCV =
% 
%    0.265951004850587
% 
% 
% miiixCV =
% 
%   -0.247886205972993
% 
%  ts: 0.0038158 r: 9.9998e-05 ei: 21170 ti: 41767
% Elapsed time is 367.339920 seconds.
% 
% maaaxJDdiff1 =
% 
%    0.056208976102692
% 
% 
% maaaxCVdiff1 =
% 
%    0.031147014705049
% 
% 
% maaaxdiff1 =
% 
%    0.079366667781697
% 
% ts: 1.2655 r: 9.9842e-05 ei: 333 ti: 632
% Elapsed time is 6.382919 seconds.
% 
% maaaxJDdiff2 =
% 
%    0.056208976102692
% 
% 
% maaaxCVdiff2 =
% 
%    0.030058296183459
% 
% 
% maaaxdiff2 =
% 
%    0.076891598822885








