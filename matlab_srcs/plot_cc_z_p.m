% plot cc & z, p value relationships

% plot cc & z 
cc=[-200:200]/200;
figure
plot(cc,fisher)
title('Correlation & Fisher-z')
xlabel({'correlation', '[-1,1]'})
ylabel('fisher-z param')
%plot cc & p
p=[];
c=[];
for freedom_idx=1:20
	dev = 1/sqrt(401/freedom_idx -3);
	fisher = 0.5.*log((1+cc)./(1-cc));
	p_tmp = abs(1-2*(1-normcdf(fisher/dev)));
	c=[c cc'];
	p=[p p_tmp'];
end	

figure

c=[c cc']
thres=ones(401,1)*0.95
p=[p thres]
% legend(threshold, 'p=0.95')
freedom_graph=plot(c, p)
title('Correlation & p value')
xlabel({'correlation','[-1,1]'})
ylabel({'p value', '[0,1]'})
