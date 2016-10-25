function auc = roc_curve(deci,label_y)
	[~,ind] = sort(deci,'descend');
	roc_y = label_y(ind);
	stack_x = cumsum(roc_y == 0)/sum(roc_y == 0);
	stack_y = cumsum(roc_y == 1)/sum(roc_y == 1);
	auc = sum((stack_x(2:length(roc_y),1)-stack_x(1:length(roc_y)-1,1)).*stack_y(2:length(roc_y),1));
        %Comment the above lines if using perfcurve of statistics toolbox
        %[stack_x,stack_y,thre,auc]=perfcurve(label_y,deci,1);
end

