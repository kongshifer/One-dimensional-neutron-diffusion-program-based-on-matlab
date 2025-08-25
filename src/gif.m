function [] = gif(DPLOT,phi_n,i_max)

num = i_max+1;% 参数a的个数，及组成GIF图的总图片个数
a_list = linspace(0, i_max, i_max+1);% 设置参数a的取值范围
x = DPLOT;% 设置自变量x的范围

pic_num = 1;
for i = 1:num  
    a = a_list(i);
    y = phi_n(:,i);
    figure(1);
    set(figure(1), 'Color', 'white');% 设置图片窗口背景颜色为白色
    plot(x, y, 'LineWidth', 2, 'Color', [0.0118, 0.0359, 0.4824], 'DisplayName', 'Number of iterations='+string(roundn(a, -2))); 
    % 绘制x,y曲线，并设置线宽，曲线颜色，曲线图例名称

    grid on;% 为绘图窗口加上网格
    set(gca,'FontSize',12,'FontName','Bookman Old Style');% 设置图片中字体的大小，样式
    xlabel('Position(cm)', 'FontSize',14);% 设置x轴标签
    xlim([min(x), max(x)]);% 设置y轴标签
    ylim([0, 1.05*max(phi_n(:,i_max+1))]);% 设置y轴显示范围
    ylabel('Neutron flux(n·cm^{-2}·s^{-1})', 'FontSize',14);% 设置y轴标签
    legend('FontSize',14, 'box', 'off');% 为图片加上图例
    title('GIF: One-dimensional problem iterative process', 'FontSize',14);% 增加图片的标题
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [A, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(A, map, 'one.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 1);
        % 将第一张图片写入one.gif文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(A, map,'one.gif','gif','WriteMode','append','DelayTime',0.8);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
end

end