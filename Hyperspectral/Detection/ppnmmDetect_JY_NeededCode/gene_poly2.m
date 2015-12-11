function [y] = gene_poly2(MPlus,alpha,b)
y=(MPlus*alpha)+ b*(MPlus*alpha).^2;
