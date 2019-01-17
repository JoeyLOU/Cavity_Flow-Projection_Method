% running script
%LDCF_ND(filename,flag,node,incre,N,Re,dt,tspan)
clc; clear; close all;

incre = 50;

% 69: 1:222
% 79: 1:186
filename = '49grid/49grid';
%LDCF_ND(filename,1,1727,incre,49,15000,1e-3,500)

node_start = 100; node_final = 1700;
plotincre = 100;
LDCF_plot(filename,node_start,node_final,plotincre,incre)
