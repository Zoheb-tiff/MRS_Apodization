%% For n peaks
tic
clearvars -except ppm BLC BL start_ppm stop_ppm

prompt = "Enter number of Metabolites: ";
n = input(prompt);

x = linspace(start_ppm, stop_ppm, length(BLC));
y = BLC;

figure, plot(ppm,BLC), title('Select the estimated Peaks')
ax = gca;
ax.XDir = 'reverse';
[estMean,estAmp] = ginput(n);
close all

figure, plot(ppm,BLC), title('Select the estimated FWHMs')
ax = gca;
ax.XDir = 'reverse';
[x1,y1] = ginput(2*n);
close all
for i = 1:n
    estFWHM(i) = abs(x1(2*i)-x1(2*i-1));
end

fun = @(X)(sum(X(:,1).*gaussian(x,X(:,2),X(:,3))) - y');
%fun = @(X)(((X(:,1).*gaussian(x,X(:,2),X(:,3)))));
lb = [];
ub = [];
x0 = [];
for i = 1:n
    lb = [lb; 0.2*estAmp(i), 0.999*estMean(i), 0.95*estFWHM(i)];
    ub = [ub; 1.01*estAmp(i), 1.001*estMean(i), 1.05*estFWHM(i)];
    x0 = [x0; estAmp(i), estMean(i), estFWHM(i)];
end% Initial Guess of the parameters


options = optimset('MaxFunEvals', 50000, 'MaxIter', 50000, 'Algorithm','trust-region-reflective', 'Display', 'off','TolX',1e-1);
[X,~,residual,~,~,~,J] = lsqnonlin(fun,x0,lb,ub,options);

Objective_fun = residual;

ans = [];
temp = 0;
for i = 1:n
    ans = [ans;X(i,1).*gaussian(x, X(i,2),X(i,3))]; % For Individual gaussians
    a = (X(i,1).*gaussian(x, X(i,2),X(i,3)));
    temp = temp + a; % For summed Gaussian
    plot(ppm, a)
    ax = gca;
    ax.XDir = 'reverse';
    lgd{i} = strcat('Gaussian',num2str(i)) ;
    hold on
end

residual = BLC' - temp;

GMM = ans;
plot(ppm,BLC)
ax = gca;
ax.XDir = 'reverse';
legend(lgd)




hold on
plot(ppm,temp,'--k','LineWidth',2)
ax = gca;
ax.XDir = 'reverse';

plot(ppm,residual)
ax = gca;
ax.XDir = 'reverse';



toc

%% CRLB calculation

J = full(J);
sigma_sq = var(residual);
FIM = (1/sigma_sq)*(real(eye(size(J,2))*(J')*(J)*eye(size(J,2))));
CRLB = diag(real(sqrt(pinv(FIM))))';
CRB = reshape(diag(real(sqrt(pinv(FIM))))',[3,n])';

%% Normalized Residual for each fit(This needs to be metabolite wise)

resid_range = zeros(n,2);

figure, plot(ppm,GMM), title('Select the tails of the gaussians')
ax = gca;
ax.XDir = 'reverse';
[estppm_resid, rr] = ginput(2*n);
resid_range = reshape(estppm_resid, [n,2])

for i = 1:size(GMM,1)
    std_peak(i) = X(i,3)/2.35;
    clearvars resid_range_idx
    temp_range_l = ppm(ppm > resid_range(i,1));
    %resid_range_idx(i,:) = find(ppm(ppm> resid_range(i,1)))
    idx_l(i) = find(ppm == temp_range_l(1));
    %clearvars resid_range_idx
    %resid_range_idx(i,:) = find(ppm(ppm< resid_range(i,2)))
    temp_range_u = ppm(ppm < resid_range(i,2));
    idx_u(i) = find(ppm == temp_range_u(end));
    clearvars resid_range_idx  
    resid = residual(idx_u(i):idx_l(i));
    fiterror(i) = (100*std(resid(:)))/(max(GMM(i,:)));
    clearvars resid
end

estFWHM = estFWHM';
fiterror = fiterror';

% resid = BLC' - temp;

% crlb = 2*(var(resid))^2/196
% area = trapz(ppm,temp)
% % relative CRLB expressed as a percentage of the presently estimated area or concentration value
% Per_CRLB = (crlb/area)*100


% D=full(J);
% F=(1/var(resid))*real(D'*D);%mean or var?
% se=real(sqrt(pinv(F)));
% plot(x,y,'r-',x,fun(X)+y,'b-')
% legend('BLC','Fitted function')


Headers = {'est_Amp','est_FWHM','est_Mean','fiterror','Calculated_FWHM'}
Final_Fit_Params = cell2table(cell(0,5),'VariableNames',Headers)
Final_Fit_Params = table(estAmp,estFWHM,estMean,fiterror,Calculated_FWHM)

writetable(Final_Fit_Params,'MRS_Fit_Params.xls')
