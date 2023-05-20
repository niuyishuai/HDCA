function setupfig(xtitle,ytitle)
% 设置坐标轴
axes1 = gca; %axes('Parent',h);
box(axes1,'on');
set(axes1,'looseInset',[0 0 0 0]);
fontsz = 13;

% 设置其余坐标区属性
set(axes1,'FontWeight','bold','FontSize',fontsz,'LineWidth',1);

% labels
xlabel(xtitle,'FontWeight','bold','FontSize',fontsz);
%ylabel(ytitle,'Interpreter','latex','FontWeight','bold','FontSize',18);
ylabel(ytitle,'FontWeight','bold','FontSize',fontsz);

%grid on;
h=gcf;
% 设置 figure属性
set(h,'InvertHardcopy','off','PaperUnits','points',...
    'Color',[1 1 1],...
    'Renderer','painters',...
    'position',[100 300 800 450]);
end