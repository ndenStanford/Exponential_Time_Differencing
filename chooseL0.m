clc
clear all
close all
figure
load('step0')
imagesc(L0)
title('initial pressure')
colorbar

figure
load('step0_mid')
imagesc(L0)
title('average pressure')
colorbar

figure
load('step0_low')
imagesc(L0)
title('bottom hole pressure')
colorbar
