
addpath(genpath('/Users/bramble/My Drive/_JPL_AMD/'))

% load spectra
A = importdata('basic_endmember_library.csv');
wavelengths = A.data(1,:);
spectra = A.data(2:end,:);
% plot spectra
plot(wavelengths,spectra)