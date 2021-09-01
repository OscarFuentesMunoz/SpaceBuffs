clear
close all
clc

a = importdata('C:\Users\super\Documents\MATLAB\Andrea Milani\deb_train\eledebtrain001.dat');
b = importdata('C:\Users\super\Documents\MATLAB\Andrea Milani\labels_train.dat');

id = 1; % change id to the corresponding debris training file number
true_sat = b(id,1);
true_conj_time = b(id,3);
true_cram = b(id,2);

debris = a(2,:);
cram_l = true_cram-10^-10;
cram_u = true_cram+10^-10;
type = 'moid';
% type = 'oe';

[cram_nom,sat,out] = debris_mother(cram_l,cram_u,debris,type);