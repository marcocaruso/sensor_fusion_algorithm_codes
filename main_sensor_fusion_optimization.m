% Main code to test the optimization loop of each SFA
% The absolute orientation errors are obtained for each unit for all the combination of parameter given as input.
% Filter codes are available as open functions in the "SFA" folder.

% When using these codes, please cite: M. Caruso, A. M. Sabatini, D. Laidig, T. Seel, M. Knaflitz, U. Della Croce, and A. Cereatti,
% “Analysis of the Accuracy of Ten Algorithms for Orientation Estimation Using Inertial and Magnetic Sensing under Optimal Conditions: One Size Does Not Fit All,” Sensors, vol. 21, no. 7, p. 2543, Apr. 2021.
% https://www.mdpi.com/1424-8220/21/7/2543

% ------------------------------------------------------------------------
% Included filters:

% GUO: S. Guo, J. Wu, Z. Wang, and J. Qian, “Novel MARG-Sensor Orientation Estimation Algorithm Using Fast Kalman Filter,” J. Sensors, vol. 2017, 2017, doi: 10.1155/2017/8542153
% LIG: G. Ligorio and A. M. Sabatini, “A novel kalman filter for human motion tracking with an inertial-based dynamic inclinometer,” IEEE Trans. Biomed. Eng., vol. 62, no. 8, pp. 2033–2043, 2015, doi: 10.1109/TBME.2015.2411431
% MAD: S. O. H. Madgwick, A. J. L. Harrison, and R. Vaidyanathan, “Estimation of IMU and MARG orientation using a gradient descent algorithm,” IEEE Int. Conf. Rehabil. Robot., vol. 2011, 2011, doi: 10.1109/ICORR.2011.5975346
% MAH: R. Mahony, T. Hamel, and J. M. Pflimlin, “Nonlinear complementary filters on the special orthogonal group,” IEEE Trans. Automat. Contr., vol. 53, no. 5, pp. 1203–1218, 2008, doi: 10.1109/TAC.2008.923738.
% MCF: MathWorks (R) Kalman filter, https://it.mathworks.com/help/releases/R2020a/fusion/ref/complementaryfilter-system-object.html?searchHighlight=complementaryfilter&s_tid=doc_srchtitle
% MKF: MathWorks (R) Complementary filter, https://it.mathworks.com/help/releases/R2020a/fusion/ref/ahrsfilter-system-object.html?searchHighlight=ahrsfilter&s_tid=doc_srchtitle
% SAB: A. M. Sabatini and S. Member, “Quaternion-Based Extended Kalman Filter for Determining Orientation by Inertial and Magnetic Sensing,” IEEE Trans. Biomed. Eng., vol. 53, no. 7, pp. 1346–1356, 2006.
% SEL: T. Seel and S. Ruppin, “Eliminating the Effect of Magnetic Disturbances on the Inclination Estimates of Inertial Sensors,” IFAC-PapersOnLine, vol. 50, no. 1, pp. 8798–8803, 2017, doi: 10.1016/j.ifacol.2017.08.1534.
% VAC: R. G. Valenti, I. Dryanovski, and J. Xiao, “Keeping a good attitude: A quaternion-based orientation filter for IMUs and MARGs,” Sensors (Switzerland), vol. 15, no. 8, pp. 19302–19330, 2015, doi: 10.3390/s150819302.
% VAK: R. G. Valenti, I. Dryanovski, and J. Xiao, “A linear Kalman filter for MARG orientation estimation using the algebraic quaternion algorithm,” IEEE Trans. Instrum. Meas., vol. 65, no. 2, pp. 467–481, 2016, doi: 10.1109/TIM.2015.2498998.

% GUO: implementation by Marco Caruso from the original implementaiton by Jin Wy available at https://github.com/zarathustr/FKF
% LIG: implementation by Angelo Maria Sabatini and Gabriele Ligorio following the original article
% MAD: implementation by Madgwick, available at https://x-io.co.uk/open-source-imu-and-ahrs-algorithms/ last accessed September, 9th 2019
% MAH: implementation by Madgwick, available at https://x-io.co.uk/open-source-imu-and-ahrs-algorithms/ last accessed September, 9th 2019
% MCF: implementation by MathWorks https://it.mathworks.com/help/releases/R2020a/fusion/ref/complementaryfilter-system-object.html?searchHighlight=complementaryfilter&s_tid=doc_srchtitle
% MKF: implementation by MathWorks https://it.mathworks.com/help/releases/R2020a/fusion/ref/ahrsfilter-system-object.html?searchHighlight=ahrsfilter&s_tid=doc_srchtitle
% SAB: implementation by Marco Caruso following the original article
% SEL: implementation by Daniel Laidig and Thomas Seel following the original article
% VAC: implementation by Marco Caruso following the original article and the C++ code avilable at https://github.com/ccny-ros-pkg/imu_tools/tree/indigo/imu_complementary_filter last accessed April 10th, 2018
% VAK: implementation by Marco Caruso following the original article
 
% ------------------------------------------------------------------------

% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

% ------------------------------------------------------------------------
%% 
clearvars
close all
clc

%% Load data
% When using this dataset, please cite:
% M. Caruso, A. M. Sabatini, M. Knaflitz, M. Gazzoni, U. D. Croce and A. Cereatti, "Orientation Estimation Through Magneto-Inertial Sensor Fusion: A Heuristic Approach for Suboptimal Parameters Tuning," in IEEE Sensors Journal, vol. 21, no. 3, pp. 3408-3419, 1 Feb.1, 2021, doi: 10.1109/JSEN.2020.3024806.
% https://ieeexplore.ieee.org/document/9201115

% dataset can be found at:
% https://ieee-dataport.org/documents/mimuopticalsassaridataset
% or at
% https://github.com/marcocaruso/mimu_optical_dataset_caruso_sassari

load('C:\Users\marco\Dropbox\caruso marco\9Febbraio\Segnali\slow_v4')
addpath(genpath('function_utilities'))
addpath(genpath('SFA'))
addpath('Optimization codes')



%% MIMU
s1 = XS1; s2 = XS2;

% Bias subtraction for each unit
s1(:,5:7) = s1(:,5:7) - mean(s1(1:5000,5:7));
s2(:,5:7) = s2(:,5:7) - mean(s2(1:5000,5:7));

% Extract the indeces correspondent to the dynamic part of the recording
indMov=[indz; indx; indy; indarb];

fs=100; % (Hz) sampling frequency

% Create time vector
t = (0:1/fs:(size(s1,1)-1)/fs)';

s1(:,1) = t;
s2(:,1) = t;
%% Stereophotogrammetry
% Refer the ground-truth orientatin to the first frame
q0  = quatconj(repmat(Qs(1,:),size(Qs,1),1));
Qs_ = quatmultiply(q0,Qs);

Qs_ = correctQuat(Qs_);

%% Init orientation quaternion
% algebraic quaternion to initialize each orientation
qin1  = quatconj(initQuaternion(-s1(1,2:4), s1(1,8:10)));
qin2  = quatconj(initQuaternion(-s2(1,2:4), s2(1,8:10)));

%% MAD
beta_ = [0 0.5 1]; % (rad/s)

[errMAD_1, errMAD_2] = EXPLORE_MAD(s1, s2, beta_, Qs_, qin1, qin2, indMov, fs);

%% MAH
kp_ = [0.1 0.5 1];   % (rad/s)
ki_ = [0.1 0.3 0.5]; % (rad/s)

[errMAH_1, errMAH_2] = EXPLORE_MAH(s1, s2, kp_, ki_, Qs_, qin1, qin2, indMov, fs);

%% GUO
stdGyr_=   [0 0.1 0.3]; % (rad/s)
stdAcc_ =  0.01;        % (m/s^2)
stdMag_ =  0.01;        % (uT)

[errGUO_1, errGUO_2] = EXPLORE_GUO(s1, s2, stdGyr_, stdAcc_, stdMag_, Qs_, indMov, qin1, qin2, fs);

%% SAB
stdGyr_ = [1e-19 0.01 0.02];    % (rad/s)
stdAcc_ = 10*9.81/1000;         % (mg)
stdMag  = 1e-3;                 % (a.u.)
thAcc_  = [0 10 200]/1000*9.81; % (mg)
thMag   = 50e-3;                % (uT)
stdBiasMag = 1e-4;              % (a.u.)

[errSAB_1, errSAB_2] = EXPLORE_SAB(s1, s2, stdGyr_, stdAcc_, stdMag, thAcc_, thMag, stdBiasMag, Qs_, qin1, qin2, indMov, fs);

%% VAC
gainAcc_  = 0.01;           % (a.u.)
gainMag_  = [0 0.01 0.5];   % (a.u.)
th1a_     = 0.05;           % (a.u.)
th2a_     = [0.05: 0.5 1];  % (a.u.)
biasalpha = 0.01;           % (a.u.)

[errVAC_1, errVAC_2] = EXPLORE_VAC(s1, s2, gainAcc_, gainMag_, th1a_, th2a_, biasalpha, Qs_, indMov, fs);

%% VAK
stdGyr_= [0 0.004];  % (rad/s)
stdAcc_ = [0 0.014]; % (m/s^2)
stdMag = 0.001;      % (uT)

[errVAK_1, errVAK_2] = EXPLORE_VAK(s1, s2, stdGyr_, stdAcc_, stdMag, Qs_, indMov, fs);

%% MCF
accGain_ = 0.01;         % (a.u.)
magGain_ =  [0 0.01 1];  % (a.u.)

[errMCF_1, errMCF_2] = EXPLORE_MCF(s1, s2, accGain_, magGain_, Qs_, qin1, qin2, indMov, fs);

%% MKF
varGyr_ =  [1e-19 0.01 0.05];  % (rad/s)^2

[errMKF_1, errMKF_2] = EXPLORE_MKF(s1, s2, varGyr_, Qs_, qin1, qin2, indMov, fs);

%% LIG

nsamp = 5000;
Ts = 1/fs;

ca_ = 1/Ts;                        % (a.u.)
cb_ = [1e-19 2]/sqrt(Ts);          % (a.u.)

stdGyr_ = deg2rad([1e-19 0.5]);    % (rad/s) 
stdAcc = 5*(9.81/1000);            % (m/s^2)      
stdMag = 1;                        % (uT) 

cm = 1;                            % magnetometer disturbance     - good for magnetically clean environments


% Set constants
% set the gravity reference frame
g  = [0; 0; 9.81];
gn = [0; 0; 1];

% set the magnetic field reference frame
localMagField = 45;

h1   = mean(s1(1:nsamp,8:10));
hn1  = h1./norm(h1);

h2 = mean(s2(1:nsamp,8:10));
hn2  = h2./norm(h2);


acc_in_1     = mean(s1(1:nsamp,2:4));
mag_in_1     = hn1*localMagField;
init1.g      = g;
init1.gn     = gn;
init1.h      = hn1*localMagField;
init1.hn     = hn1;
init1.normh  = norm(h1);
init1.acc_in = acc_in_1;
init1.mag_in = mag_in_1;


acc_in_2     = mean(s2(1:nsamp,2:4));
mag_in_2     = hn2*localMagField;
init2.g      = g;
init2.gn     = gn;
init2.h      = hn2*localMagField;
init2.normh  = norm(h2);
init2.hn     = hn2;
init2.acc_in = acc_in_2;
init2.mag_in = mag_in_2;

[errLIG_1,errLIG_2] = EXPLORE_LIG(s1, s2, stdGyr_, stdAcc, stdMag, ca_, cb_, cm, Qs_, init1, init2, indMov, localMagField, fs);

%% SEL
tauAcc_= [1e-19 1];
tauMag_= [1e-19 3];
zeta_= 0;
accRating_= 1;

[errSEL_1,errSEL_2] =  EXPLORE_SEL(s1,s2,tauAcc_, tauMag_, zeta_, accRating_,Qs_,qin1,qin2,indMov,fs);