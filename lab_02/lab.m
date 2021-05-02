function main()
    X=[-5.05,-5.74,-6.39,-5.01,-4.94,-6.32,-4.73,-5.44,-4.79,-5.40,...
       -4.50,-3.43,-5.21,-5.22,-5.07,-5.51,-4.45,-5.24,-6.50,-4.99,...
       -5.42,-3.30,-5.06,-5.38,-6.48,-3.31,-5.56,-6.50,-5.22,-5.68,...
       -4.80,-6.67,-4.95,-5.32,-5.68,-6.32,-3.72,-4.59,-6.33,-5.03,...
       -4.49,-4.80,-6.04,-6.21,-3.60,-3.93,-5.89,-5.29,-7.41,-3.73,...
       -6.61,-4.34,-5.99,-5.24,-4.08,-4.68,-5.38,-6.38,-4.66,-3.67,...
       -4.61,-4.54,-4.51,-5.43,-6.47,-5.31,-4.30,-6.32,-5.82,-3.44,...
       -5.92,-4.76,-4.45,-3.52,-4.91,-5.65,-5.02,-5.00,-5.26,-4.98,...
       -6.16,-6.21,-4.42,-6.20,-5.84,-5.58,-5.34,-5.21,-5.78,-7.80,...
       -5.21,-4.79,-4.53,-4.78,-6.39,-7.04,-4.82,-5.53,-3.52,-6.24,...
       -3.58,-5.01,-5.79,-4.80,-6.04,-5.15,-7.03,-4.71,-4.38,-5.77,...
       -4.05,-5.76,-5.86,-6.45,-4.81,-5.68,-7.48,-3.97,-5.16,-3.48];
   
    % Уровень доверия
    gamma = 0.9;
    % Объем выборки 
    n = length(X);
    % Точечная оценка матожидания
    mu = mean(X);
    % Точечная оценка дисперсии
    s2 = var(X);
    
    % Нижняя граница доверительного интервала для матожидания
    muLow = findMuLow(n, mu, s2, gamma);
    % Верхняя граница доверительного интервала для матожидания
    muHigh = findMuHigh(n, mu, s2, gamma);
    % Нижняя граница доверительного интервала для дисперсии
    s2Low = findS2Low(n, s2, gamma);
    % Верхняя граница доверительного интервала для дисперсии
    s2High = findS2High(n, s2, gamma);
    
    % Вывод полученных ранее значений
    fprintf('mu = %.3f\n', mu);
    fprintf('S2 = %.3f\n', s2);
    fprintf('muLow = %.3f\n', muLow);
    fprintf('muHigh = %.3f\n', muHigh);
    fprintf('s2Low = %.3f\n', s2Low);
    fprintf('s2High = %.3f\n', s2High);
    
    % Создание массивов точченых оценок
    muArray = zeros(1, n);
    s2Array = zeros(1, n);
    % Создание массивов границ доверительных интервалов
    muLowArray = zeros(1, n);
    muHighArray = zeros(1, n);
    s2LowArray = zeros(1, n);
    s2HighArray = zeros(1, n);
    
    % Цикл от 1 до n
    for i = 1 : n
        mu = mean(X(1:i));
        s2 = var(X(1:i));
        % Точечная оценка матожидания
        muArray(i) = mu;
        % Точечная оценка дисперсии
        s2Array(i) = s2;
        % Нижняя граница доверительного интервала для матожидания
        muLowArray(i) = findMuLow(i, mu, s2, gamma);
        % Верхняя граница доверительного интервала для матожидания
        muHighArray(i) = findMuHigh(i, mu, s2, gamma);
        % Нижняя граница доверительного интервала для дисперсии
        s2LowArray(i) = findS2Low(i, s2, gamma);
        % Верхняя граница доверительного интервала для дисперсии
        s2HighArray(i) = findS2High(i, s2, gamma);
    end
    
    % Построение графиков
    plot(1 : n, [(zeros(1, n) + mu)', muArray', muLowArray', muHighArray']);
    xlabel('n');
    ylabel('y');
    legend('$\hat \mu(\vec x_N)$', '$\hat \mu(\vec x_n)$', ...
        '$\underline{\mu}(\vec x_n)$', '$\overline{\mu}(\vec x_n)$', ...
        'Interpreter', 'latex', 'FontSize', 18);
    figure;
    plot(1 : n, [(zeros(1, n) + s2)', s2Array', s2LowArray', s2HighArray']);
    xlabel('n');
    ylabel('z');
    legend('$\hat S^2(\vec x_N)$', '$\hat S^2(\vec x_n)$', ...
        '$\underline{\sigma}^2(\vec x_n)$', '$\overline{\sigma}^2(\vec x_n)$', ...
        'Interpreter', 'latex', 'FontSize', 18);
end

% Функция поиска нижней границы доверительного интервала для матожидания
function muLow = findMuLow(n, mu, s2, gamma)
    muLow = mu - sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

% Функция поиска верхней границы доверительного интервала для матожидания
function muHigh = findMuHigh(n, mu, s2, gamma)
    muHigh = mu + sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

% Функция поиска нижней границы доверительного интервала для дисперсии
function s2Low = findS2Low(n, s2, gamma)
    s2Low = ((n - 1) * s2) / chi2inv((1 + gamma) / 2, n - 1);
end

% Функция поиска верхней границы доверительного интервала для дисперсии
function s2High = findS2High(n, s2, gamma)
    s2High = ((n - 1) * s2) / chi2inv((1 - gamma) / 2, n - 1);
end