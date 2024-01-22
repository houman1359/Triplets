




% Plot 2D margins
figure('Name','2D margins of mixed vine copula PDF','Position',[0,0,1600,1000]);
subplot(2,3,1);
margin12 = reshape(sum(sum(p,4),3),[length(x1gv),length(x2gv)]);
imagesc(x1gv,x2gv,margin12');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 2');
subplot(2,3,2);
margin13 = reshape(sum(sum(p,4),2),[length(x1gv),length(x3gv)]);
imagesc(x1gv,x3gv,margin13');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 3');
subplot(2,3,3);
margin14 = reshape(sum(sum(p,3),2),[length(x1gv),length(x4gv)]);
imagesc(x1gv,x4gv,margin14');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 4');
subplot(2,3,4);
margin23 = reshape(sum(sum(p,4),1),[length(x2gv),length(x3gv)]);
imagesc(x2gv,x3gv,margin23');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 2');
ylabel('Margin 3');
subplot(2,3,5);
margin24 = reshape(sum(sum(p,3),1),[length(x2gv),length(x4gv)]);
imagesc(x2gv,x4gv,margin24');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 2');
ylabel('Margin 4');
subplot(2,3,6);
margin34 = reshape(sum(sum(p,2),1),[length(x3gv),length(x4gv)]);
imagesc(x3gv,x4gv,margin34');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 3');
ylabel('Margin 4');
fprintf('\n');





%% Test sampling
disp('Sampling from mixed copula vine...');
% Draw samples
cases = 1000;
x = mixedvinernd(vineest,cases);
% Plot samples in 2D
figure('Name','Mixed vine copula samples in 2D','Position',[0,0,1600,1000]);
subplot(2,3,1);
scatter(x(:,1),x(:,2),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
subplot(2,3,2);
scatter(x(:,1),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 3');
subplot(2,3,3);
scatter(x(:,1),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 4');
subplot(2,3,4);
scatter(x(:,2),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 3');
subplot(2,3,5);
scatter(x(:,2),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 4');
subplot(2,3,6);
scatter(x(:,3),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 3');
ylabel('Margin 4');
% Plot samples in 3D
figure('Name','Mixed vine copula samples in 3D','Position',[0,0,1000,1000]);
subplot(2,2,1);
scatter3(x(:,1),x(:,2),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
zlabel('Margin 3');
subplot(2,2,2);
scatter3(x(:,1),x(:,2),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
zlabel('Margin 4');
subplot(2,2,3);
scatter3(x(:,1),x(:,3),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 3');
zlabel('Margin 4');
subplot(2,2,4);
scatter3(x(:,2),x(:,3),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 3');
zlabel('Margin 4');
fprintf('\n');




