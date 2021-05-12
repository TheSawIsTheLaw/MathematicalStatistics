function main()
    X=[-14.34,-16.97,-14.09,-14.74,-16.69,-13.85,-15.55,-14.62,-13.30,-15.52,...
       -14.75,-16.51,-17.15,-16.87,-15.06,-13.60,-14.48,-14.71,-14.17,-13.88,...
       -14.55,-15.37,-14.81,-16.05,-17.06,-15.86,-15.12,-15.98,-14.16,-15.81,...
       -15.06,-16.19,-16.22,-16.19,-14.87,-15.62,-15.86,-15.25,-16.34,-14.44,...
       -14.72,-15.17,-15.24,-14.44,-15.93,-14.87,-16.53,-15.76,-15.12,-12.91,...
       -16.06,-16.06,-14.89,-15.57,-13.59,-16.84,-13.88,-14.33,-15.45,-16.58,...
       -16.05,-14.34,-13.55,-16.78,-14.15,-14.28,-14.40,-13.98,-16.23,-15.35,...
       -14.77,-15.61,-15.59,-15.64,-14.76,-17.18,-15.13,-15.01,-14.21,-13.91,...
       -16.55,-15.44,-14.03,-16.44,-15.57,-15.07,-16.28,-16.30,-15.74,-14.03,...
       -14.85,-15.73,-15.81,-14.42,-14.14,-15.14,-15.49,-16.42,-14.22,-14.20,...
       -17.17,-15.82,-14.96,-14.75,-14.98,-13.64,-14.00,-17.29,-14.51,-16.18,...
       -15.70,-15.07,-14.28,-14.55,-13.85,-15.36,-15.74,-14.61,-16.32,-15.34];
    
    % 1)
    % а) Максимальное и минимальное значения
    Mmax = max(X);
    Mmin = min(X);
    fprintf("\nа) Mmax (максимальное значение) = %f; Mmin (минимальное значенение) = %f", Mmax, Mmin);
    
    % б) Размах
    R = Mmax - Mmin;
    fprintf("\nб) R (размах) = %f", R);
    
    % в) Оценки
    mu = sum(X) / length(X);
    s2 = sum((X - mean(X)).^2) / (length(X) - 1);
    fprintf("\nв)mu (оценка математического ожидания) = %f; s^2 (оценка дисперсии) = %f", mu, s2);
    
    % г) Группировка значений выборки
    % Нахождение количества интервалов
    m = floor(log2(length(X))) + 2;
    fprintf("\nг)Группировка значений выборки в m = [log2 n] + 2 интервала: m = %f\n", m);
    
    % Разбиение выборки на интервалы от min до max, с помощью BinLimits
    % объединяем только те значения, которые находятся в интервале от
    % минимума до максимума
    [counts, edges] = histcounts(X, m, 'BinLimits', [min(X), max(X)]);
    
    for i = 1: length(counts)
        fprintf("[%f : %f] - %d\n", edges(i), edges(i + 1), counts(i));
    end
    
    % д) Построение гистограмы
    hist = histogram();
    hist.BinEdges = edges;
    hist.BinCounts = counts / length(X) / ((max(X) - min(X)) / m);
    
    hold on; % Продолжаем работать с той же системой
    
    % График функции плотности рапределения вероятностей нормальной случайной величины
    delta = R/m;
    sigma = sqrt(s2);
    Xn = min(X):delta/20:max(X);
    Y = normpdf(Xn, mu, sigma);
    plot(Xn, Y, 'blue');
    
    % e)
    figure;
    [yy, xx] = ecdf(X);
    stairs(xx, yy);
    
%     uniques = unique(X);
%     count = histcounts(X, uniques);
%     for i = 2 : (length(count))
%         count(i) = count(i) + count(i - 1);
%     end
%     count = [0 count];
%     stairs(uniques, count / length(X));
    
    hold on;
    
    delta = R/m;
    Xn = min(X):delta/20:max(X);
    Y = normcdf(Xn, mu, s2);
    plot(Xn, Y, 'black');
    
end
