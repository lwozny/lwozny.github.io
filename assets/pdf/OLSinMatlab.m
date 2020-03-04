clear; % czyszczenie pami?ci 

importeddata = csvread('TableF2-2.csv'); % import danych z pliku csv w OCTAVE

data = importeddata(2:end,:); % stworzenie macierzy data z pominieciem pierwszego wiersza macierzy importeddata 

y = log(data(:,2)./data(:,3)); % stworzenie zmiennej objasnianej y poprzez policzenie logarytmu z ilorazu 2giej i 3ciej kolumny macierzy data. 

G = log(data(:,4:7)); % wybranie logarytmów zmiennych objasnianych z kolumn od 4 do 7 

z = ones(52,1); % stworzenie wektora jedynek (do policzenia sta?ej w regresji)

X = [z , G]; % stworzenie wektora zmiennych obja?nianych poprzez po??czenie z oraz G

inv(X'*X)*X'*y % wyznaczanie parametrów regresji (elastyczno?ci) metod? OLS (z definicji)

%alternatywnie mozemy skorzystac z funkcji regress z pakietu statistics.
%?adowanie pakietu w OCTAVE: pkg load statistics

[b, bcof, r, rconf, stat] = regress(y,X);

b

stat

plot(rconf);