function [quat, bias, error] = runOriEstIMU(acc, gyr, mag, rate, tauAcc, tauMag, zeta, accRating, qin)
%RUNORIESTIMU Convenience function to run OriEstIMU
%   acc, gyr, mag: input data as Nx3 matrices, for 6D sensor fusion use mag=zeros(size(acc))
%   rate: sampling frequeny in Hz
%   tauAcc, tauMag, zeta, accRating: tuning parameters
%   qin: initial orientation quaternion

    N = size(acc, 1);
    quat = zeros(N, 4);
    bias = zeros(N, 3);
    error = zeros(N, 2);
    state = [qin zeros(1, 14)];
    for i=1:N
        [state, errorAngleIncl, errorAngleAzi] = CSG_OriEst_IMU(state, acc(i,:), gyr(i,:), mag(i,:), rate, tauAcc, tauMag, zeta, accRating);
        quat(i,:) = state(1:4);
        bias(i,:) = state(5:7);
        error(i,:) = [errorAngleIncl, errorAngleAzi];
    end
end
