%% Ajuda
% doc <ordre>
doc plot; % s'obre una finestra nova
help plot; % l'ajuda apareix en la línia d'ordres

%% Línies d'ordres
1 + 1 % un ordre
2 + 2; 4 * 5 % dos ordres separats per ;
2 + 3 + ...
    4 % un ordre escrit en més d'una línia
%% Comentaris

% En línia pròpia, abans o després dels ordres
% Calcula 1 + 1!
1 + 1
% He calculat 1 + 1.

% En la mateixa línia amb els ordres
2 * 3 % calcula 2 * 3

%% Supressió d'eixida

4 * 5 % imprimeix el resultat en la línia d'ordres
4 * 5; % no imprimeix res

%% Nombres

1 % nombre sencer
1.2 % nombre real
1e-3 % potència de 10: 10^N => 1eN, N ha de ser sencer
% 1e2.5 -> dóna error
2.5e12 % 2.5 * 10^12
1 + 3i % nombre complexe

%% Formats

format short; % per defecte
1.000000000005

format long; % recomanat
1.000000000005

format short e; % format científic curt
1.000000000005

format long e; % format científic llarg
1.000000000005

format long; % tornar a usar long

%% Impressió
fprintf('Hello world!') % imprimeix, però no passa a la seguent línia
fprintf('\n'); % \n -> comença nova línia
fprintf('What time is it?\nI do not know!\n');

fprintf('%d * %d = %d\n', 2, 3, 2 * 3); % %d és el formatador per als nombres sencers
fprintf('%d / %d = %g\n', 5, 2, 5/2); % %g es pot usar per als reals

fprintf('%.4d\n', 3); % imprimir amb zeros si calen
fprintf('%4d\n', 3); % imprimir amb espais blancs si calen
fprintf('%.12f\n', 4/3); % imprimir amb 12 xifres després de la coma
fprintf('%.12g\n', 4/3); % imprimir amb 12 xifres significants

fprintf('My name is %s.\n', 'John'); % %s es pot usar per als strings
fprintf('My name is %6s.\n', 'John'); % imprimir el nom amb espais blancs si calen

% impressió formatada
fprintf('ITEM   \t QUANTITY\n'); % \t és el tabulador
fprintf('-----------------\n');
fprintf('%6s \t %3d\n', 'Books', 12);
fprintf('%6s \t %3d\n', 'Food', 120);
fprintf('%6s \t %3d\n', 'Cars', 2);
%% Variables

x = 1 % assigna el valor 1 a la variable x i l'imprimeix
y = 1; % assigna el valor sense imprimir (recomanat)

hello = 'Hello world!'; % una variable de tipus cadena (string)

%% Operadors

x = 2;
y = 4;

% Operadors aritmètics:

x + y; x + 3; 1 + 2 % suma
x - y; -3 - 4 % resta
x * y; 2 * 3 % multiplicació
x / y; 5 / 2.8 % divisió
x ^ y; 3 ^ 3 % potenciació

% Atenció:
1 / 0 % -> Inf: divisió entre 0
0 / 0 % -> NaN: "not a number"
Inf / Inf % NaN

% Operadors relacionals:

x == y % igual
x ~= y % no igual
3 < 4 % inferior
4 > 3 % més gran
5 <= 6 % inferior o igual
7 >=6 % més gran o igual

% Parèntesi: el que està dins d'un parentesi es calcula primer
3 + 4 * 5 ~= (3 + 4) * 5

%% funcions matemàtiques
sqrt(4) % arrel quadrada
nthroot(64, 3) % arrel enèsima

pi % el valor de π
cos(pi)
sin(pi/2)
tan(pi)
cot(0) % -> Inf

acos(0)
acos(cos(pi))
asin(0)
atan(1)
acot(1)

exp(1) % e^1
log(10) % ln
log10(10) % log
%% Vectors (matemàtiques)

[1 2 3] % vector fila
[1, 2, 3] % vector fila (recomanat)

[-12; 13.45; 8] % vector columna

% trasposta:
vec_fila = [1, 2, 3, 4, 5];
vec_fila'

% producte escalar
a = [2, 3, 4];
b = [0, -1, 3];
dot(a, b) % a . b
% a * b  -> dóna error, cal multiplicar un vector fila amb un vector
% columna
a * b'

c = [10; 12; -1]
dot(a, c)
a * c

% producte vectorial
cross(a, c)

% norma
norm(a)
norm(b)
norm(c)

%% Matrius

A = [1, 2, 3;
    4, 5, 6;
    7, 8, 9];

b = [4, 5, 6];
c = [4; 5; 6];

D = [0, 2;
    -10, 12;
    1, 0.22];

% A * b -> dóna error
b * A
A * c
% c * A -> dóna error
A * D
% D * A -> dóna error

% determinant
det(A) % matriu singular
E = [1 2; 8 9];
det(E)

% inversa
inv(A) % avisa que la matriu pot ser singular
inv(E)

%% Vectors (programació)

W = [1 2 3 4 5]; % vector amb 5 elements

W(1) % primer element
W(3) % tercer element
W(end) % l'ultim element

W(4) = -10 % modificar el cuart element
W

Z = 1:30; % genera el vector 1, 2, 3, ...., 30
Z2 = 1:0.5:30; % general el vector 1, 1.5, 2, 2.5, 3, 3.5,  ... , 29.5, 30

length(Z) % longitud d'un vector
length(Z2)

Z(length(Z))
Z(length(Z)) == Z(end) % l'ultim element

% operacions vectoritzades

a = 1:10;
b = a + 3; % suma 3 a cada element
c = a - 4;
d = a * 2;

e = a .* d; % .* : multiplica els elements corresponents
% e = a * d -> dóna error
f = a ./ e; % ./ : divideix els elements corresponents

g = 1:12;

% a .* g -> dóna error, la longitud dels vectors ha de coincidir

% g^2 -> dóna error
h = g .^ 2; % potenciació de cada element

cos(h) % calcual el cosinus de cada element

x = -1:0.1:1;
x(22) = 0; % afegir més elements a un vector

%% Estructures de control

% if-then-else
a = 3;
if (a == 3) % comprova si la condició dóna 1 (veritat)
    b = a + 1; % executa aquest grup si veritat
else
    c = 0; % si no és veritat, executa aquest grup
    b = a - 1;
end

% for
aa = 0;
for i = 1:10 % executa el bloc d'ordres per als valors 1, 2, 3, ..., 10
    aa = aa + 1;
end
aa

bb = 0;
for i = -10:-1 % executa el bloc d'ordres per als valors -10, -9, ..., -1
    bb = bb + 1;
end
bb

% revers
for j = 10:-1:1
    bb = bb - 1;
end
bb

% bucles anidats
cc = 0;
for i = 1:10 % executa el bloc 10 vegades
    % bloc de bucle i
    for j = 1:3 % executa el bloc 3 vegades
        cc = cc + i * j;
        cc = cc * 2;
    end
end
cc

% while
a = 3; b = 0;
while(a >= 0)
    b = b + a;
    a = a - 1;
end
a
b

% exemple: imprimir una tabla usant el bucle for
fprintf('%3s | %3s || %3s\n', 'A', 'B', 'A+B');
fprintf('----------------\n');
A = 3; B = 8;
for i = 1:5
    fprintf('%3d | %3d || %3d\n', A, B, A + B);
    A = A + 3 * B;
    B = B - A;
end

%% Funcions

f = @(x) x + 1; % f(x) = x + 1
f(3)
a = 3; f(a)

f2 = inline('x + 1'); % altra manera de definir una funció
f2(3)

g = @(a, b, c) a + b + c; % una funció de tres variables
g(1, 2, 3)

h = @(x) [cos(x), sin(x)]; % una funció que retorna un vector
h(0)
h(pi)

myfunc(1, 4, 5) % crida la funció declarada en el fixter myfunc.m

%% Representació gràfica (ràpida)
ezplot('x^2'); % representa una funció automàticament
ezplot('x^3', [-0.5, 0.5]); % controla l'interval en x

% defineix una funció i després representa-la
f = @(x) x.^2 + sin(x);
ezplot(f);

% representa una funció implícita
ezplot('x^2 - y^2 = 2');
%% Representació gràfica (estàndard)
x = 1:10;
y = x.^2;
plot(x, y); % pintar y = x^2
plot(x, y, '.'); % pintar amb els punts
plot(x, y, 'b*'); % pintar amb els asterics blaus
plot(x, y, 'kv-'); % pintar amb els triangles negres connectats

z = x.^2 / 2;
plot(x, y, 'b*', x, z, 'kv'); % pintar dos parelles de vectors
plot(x, y, 'b*', y, z, 'r-');

w = 1:9;
% plot(x, w); -> dóna error, els vectors han de tindre la mateixa longitud

% controlar el tamany i els colors
plot(x, y, 'bo-', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'g', ...
    'MarkerFaceColor', 'r');

% mostrar els eixos i el títol
xlabel('x');
ylabel('y');
title('y = x^2');

% canviar l'eix a logarítmic
semilogy(x, y); % només l'eix y
loglog(x, y); % amdos

% múltiples gráfics
subplot(2, 2, 1); plot(x, y); title('Lloc 1');
subplot(2, 2, 2); plot(x, -y, 'r'); title('Lloc 2');
subplot(2, 2, 3); plot(-x, y, 'g'); title('Lloc 3');
subplot(2, 2, 4); plot(-x, -y, 'b'); title('Lloc 4');

% guardar en un fitxer PNG
print('-dpng', 'plot.png');

% guardar en un fitxer DF
print('-dpdf', 'plot.pdf');

%% Llegir dades d'un fitxer
gdp = readtable('GDP.csv');
summary(gdp)
plot(gdp.Year, gdp.PIB_euros, 'ko-');
xlabel('Year'); ylabel('GDP of Spain in euros');
title('GDP of Spain');

%% Guardar en un fitxer
gdp.PIB_triillions = gdp.PIB_euros / 1e12;
plot(gdp.Year, gdp.PIB_triillions, 'ko-');
xlabel('Year'); ylabel('GDP of Spain in triillions of euros');
title('GDP of Spain');
writetable(gdp, 'GDP_calc.csv');
type('GDP_calc.csv')

%% Mesura de temps bàsica
X = 1:100000; % genera un vector
tic; % inici de la mesura de temps
Y = sin(X);
toc % fin de la mesura

%% Mesura de temps avançada
num_measurements = 1000; % repetim 1000 vegades
t = zeros(1, num_measurements);
for i = 1:num_measurements
    tic;
    Y = sin(X);
    t(i) = toc;
end
fprintf("Temps d'execucio = %.4e +- %.4e segons\n", mean(t), std(t));

%% Ordres per generar vectors i matrius
v = zeros(1, 10); % genera un vector fila de 10 zeros
w = zeros(4, 1); % genera un vector columna de 4 zeros
A = zeros(5, 4); % genera una matriu de 5 files i 4 columnes

v2 = ones(1, 10); % el mateix, pero amb 1 en comptes de 0

ident = eye(10); % la matriu identitat de dimensions 10x10

% Aquests comandaments es poden combinar:

ones(10, 10) + 3 * eye(10)

%% Càlcul simbòlic

% atenció: sempre s'han de declarar les variables simbòliques:
syms x; % declara x
syms f(x); % declara f(x)
f(x) = sin(x);

diff(f(x)) % derivada
int(f(x), x) % integral
solve(f(x) == pi, x) % resoldre una equació

% resoldre una equació diferencial
syms a g(a)
dsolve(diff(g(a), a) == 1/a, a)