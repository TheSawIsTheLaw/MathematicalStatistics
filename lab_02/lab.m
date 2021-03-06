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
   
    % Уровень доверия
    gamma = 0.9;
    % Объем выборки 
    n = length(X);
    % Точечная оценка мат. ожидания
    mu = mean(X);
    % Точечная оценка дисперсии
    s2 = var(X);
    
    % Нижняя граница доверительного интервала для мат. ожидания
    muBot = findMuBot(n, mu, s2, gamma);
    % Верхняя граница доверительного интервала для мат. ожидания
    muTop = findMuTop(n, mu, s2, gamma);
    % Нижняя граница доверительного интервала для дисперсии
    s2Bot = findS2Bot(n, s2, gamma);
    % Верхняя граница доверительного интервала для дисперсии
    s2Top = findS2Top(n, s2, gamma);
    
    % Вывод полученных ранее значений
    fprintf('mu (Точечная оценка математического ожидания) = %.3f\n', mu);
    fprintf('S2 (Точечная оценка дисперсии) = %.3f\n', s2);
    fprintf('muBot (нижняя граница доверительного интервала для математического ожидания) = %.3f\n', muBot);
    fprintf('muTop (верхняя граница -//-) = %.3f\n', muTop);
    fprintf('s2Bot (нижняя граница доверительного интервала для дисперсии) = %.3f\n', s2Bot);
    fprintf('s2Top (верхняя граница -//-) = %.3f\n', s2Top);
    
    % Создание массивов точечных оценок
    muArray = zeros(1, n);
    s2Array = zeros(1, n);
    % Создание массивов границ доверительных интервалов
    muBotArray = zeros(1, n);
    muTopArray = zeros(1, n);
    s2BotArray = zeros(1, n);
    s2TopArray = zeros(1, n);
    
    for i = 1 : n
        mu = mean(X(1:i));
        s2 = var(X(1:i));
        % Точечная оценка матожидания
        muArray(i) = mu;
        % Точечная оценка дисперсии
        s2Array(i) = s2;
        % Нижняя граница доверительного интервала для матожидания
        muBotArray(i) = findMuBot(i, mu, s2, gamma);
        % Верхняя граница доверительного интервала для матожидания
        muTopArray(i) = findMuTop(i, mu, s2, gamma);
        % Нижняя граница доверительного интервала для дисперсии
        s2BotArray(i) = findS2Bot(i, s2, gamma);
        % Верхняя граница доверительного интервала для дисперсии
        s2TopArray(i) = findS2Top(i, s2, gamma);
    end
    
    % Построение графиков
    plot(1 : n, [(zeros(1, n) + mu)', muArray', muBotArray', muTopArray']);
    xlabel('n');
    ylabel('y');
    legend('$\hat \mu(\vec x_N)$', '$\hat \mu(\vec x_n)$', ...
        '$\underline{\mu}(\vec x_n)$', '$\overline{\mu}(\vec x_n)$', ...
        'Interpreter', 'latex', 'FontSize', 18);
    figure;
    plot(1 : n, [(zeros(1, n) + s2)', s2Array', s2BotArray', s2TopArray']);
    xlabel('n');
    ylabel('z');
    legend('$\hat S^2(\vec x_N)$', '$\hat S^2(\vec x_n)$', ...
        '$\underline{\sigma}^2(\vec x_n)$', '$\overline{\sigma}^2(\vec x_n)$', ...
        'Interpreter', 'latex', 'FontSize', 18);
end

% Функция поиска нижней границы доверительного интервала для матожидания
function muBot = findMuBot(n, mu, s2, gamma)
    muBot = mu - sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

% Функция поиска верхней границы доверительного интервала для матожидания
function muTop = findMuTop(n, mu, s2, gamma)
    muTop = mu + sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

% Функция поиска нижней границы доверительного интервала для дисперсии
function s2Bot = findS2Bot(n, s2, gamma)
    s2Bot = ((n - 1) * s2) / chi2inv((1 + gamma) / 2, n - 1);
end

% Функция поиска верхней границы доверительного интервала для дисперсии
function s2Top = findS2Top(n, s2, gamma)
    s2Top = ((n - 1) * s2) / chi2inv((1 - gamma) / 2, n - 1);
end
