function main()
    X=[-5.05,-5.74,-6.39,-5.01,-4.94,-6.32,-4.73,-5.44,-4.79,-5.40,-4.50,-3.43,-5.21,-5.22,-5.07,-5.51,-4.45,-5.24,-6.50,-4.99,-5.42,-3.30,-5.06,-5.38,-6.48,-3.31,-5.56,-6.50,-5.22,-5.68,-4.80,-6.67,-4.95,-5.32,-5.68,-6.32,-3.72,-4.59,-6.33,-5.03,-4.49,-4.80,-6.04,-6.21,-3.60,-3.93,-5.89,-5.29,-7.41,-3.73,-6.61,-4.34,-5.99,-5.24,-4.08,-4.68,-5.38,-6.38,-4.66,-3.67,-4.61,-4.54,-4.51,-5.43,-6.47,-5.31,-4.30,-6.32,-5.82,-3.44,-5.92,-4.76,-4.45,-3.52,-4.91,-5.65,-5.02,-5.00,-5.26,-4.98,-6.16,-6.21,-4.42,-6.20,-5.84,-5.58,-5.34,-5.21,-5.78,-7.80,-5.21,-4.79,-4.53,-4.78,-6.39,-7.04,-4.82,-5.53,-3.52,-6.24,-3.58,-5.01,-5.79,-4.80,-6.04,-5.15,-7.03,-4.71,-4.38,-5.77,-4.05,-5.76,-5.86,-6.45,-4.81,-5.68,-7.48,-3.97,-5.16,-3.48];
    
    %1)
    %а)
    Mmax = max(X);
    Mmin = min(X);
    fprintf("\nа) Mmax (максимальное значение) = %f; Mmin (минимальное значенение) = %f", Mmax, Mmin);
    
    %б)
    R = Mmax - Mmin;
    fprintf("\nб) R (размах) = %f", R);
    
    %в)
    mu = sum(X) / length(X);
    s2 = sum((X - mean(X)).^2) / (length(X) - 1);
    fprintf("\nв)mu (оценка математического ожидания) = %f; s^2 (оценка дисперсии) = %f", mu, s2);
    
    %г)
    
    %д)
    
    %е)