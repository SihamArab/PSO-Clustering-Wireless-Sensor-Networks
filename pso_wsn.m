function out = pso_wsn(varargin)
% PSO_WSN  Clustering par PSO pour reseaux de capteurs sans fil, compare a LEACH.
%
%   pso_wsn()                              parametres du rapport, 50 executions
%   pso_wsn('nRuns', 3, 'iters', 20)       surcharge de parametres, par paires nom, valeur
%
% Modele radio du premier ordre, MAC ideal, reseau statique a un saut.
% Ecrit results/summary.csv et, si l'affichage est disponible, des figures.
% Sans toolbox.
P = wsn_params();
for i = 1:2:numel(varargin)
    P.(varargin{i}) = varargin{i + 1};
end
protos = {'leach', 'pso'};
R = P.maxRounds;
for a = 1:2
    S.(protos{a}).alive = zeros(R, P.nRuns);
    S.(protos{a}).Eavg = zeros(R, P.nRuns);
    S.(protos{a}).Evar = zeros(R, P.nRuns);
    S.(protos{a}).Econs = zeros(R, P.nRuns);
    S.(protos{a}).packets = zeros(R, P.nRuns);
    S.(protos{a}).FND = nan(P.nRuns, 1);
    S.(protos{a}).LND = nan(P.nRuns, 1);
end
for run = 1:P.nRuns
    set_seed(run);
    net = init_network(P);
    for a = 1:2
        set_seed(1000 + run);
        res = run_protocol(protos{a}, P, net);
        s = S.(protos{a});
        s.alive(:, run) = res.aliveCnt;
        s.Eavg(:, run) = res.Eavg;
        s.Evar(:, run) = res.Evar;
        s.Econs(:, run) = res.Econs;
        s.packets(:, run) = res.packets;
        s.FND(run) = res.FND;
        s.LND(run) = res.LND;
        S.(protos{a}) = s;
    end
    fprintf('run %d/%d  LEACH FND=%g LND=%g | PSO FND=%g LND=%g\n', run, P.nRuns, ...
        S.leach.FND(run), S.leach.LND(run), S.pso.FND(run), S.pso.LND(run));
end

if ~exist('results', 'dir'), mkdir('results'); end
fid = fopen('results/summary.csv', 'w');
fprintf(fid, 'protocol,metric,mean,std\n');
for a = 1:2
    s = S.(protos{a});
    m = {'FND', s.FND; 'LND', s.LND; 'LND_minus_FND', s.LND - s.FND; ...
         'total_energy_consumed_J', sum(s.Econs, 1)'; ...
         'mean_packets_to_BS_per_round', (sum(s.packets, 1) ./ max(sum(s.alive > 0, 1), 1))'};
    for j = 1:size(m, 1)
        v = m{j, 2};
        fprintf(fid, '%s,%s,%.6g,%.6g\n', protos{a}, m{j, 1}, nanmean_(v), nanstd_(v));
    end
end
fclose(fid);
fprintf('Resume ecrit dans results/summary.csv\n');
out = S;

try
    fig = figure('visible', 'off');
    hold on
    plot(mean(S.leach.alive, 2), 'r'); plot(mean(S.pso.alive, 2), 'b');
    xlabel('Tour'); ylabel('Noeuds vivants'); legend('LEACH', 'PSO'); grid on
    print(fig, 'results/alive_nodes.png', '-dpng');
    close(fig)
    fig = figure('visible', 'off');
    hold on
    plot(mean(S.leach.Eavg, 2), 'r'); plot(mean(S.pso.Eavg, 2), 'b');
    xlabel('Tour'); ylabel('Energie residuelle moyenne, J'); legend('LEACH', 'PSO'); grid on
    print(fig, 'results/residual_energy.png', '-dpng');
    close(fig)
catch err
    warning('Figures non ecrites : %s', err.message);
end
end

function m = nanmean_(v), v = v(~isnan(v)); if isempty(v), m = NaN; else m = mean(v); end, end
function s = nanstd_(v), v = v(~isnan(v)); if numel(v) < 2, s = NaN; else s = std(v); end, end

function P = wsn_params()
% WSN_PARAMS  Parametres de simulation, tels que decrits dans le rapport
% (tableau 4.1 et section 4.3.4). Voir le README pour les choix qui ne sont
% pas fixes par le rapport.

% Reseau
P.area      = 100;          % cote du champ, m
P.N         = 100;          % nombre de noeuds
P.nAdv      = 20;           % noeuds avances, les autres sont normaux
P.E0normal  = 0.5;          % J
P.E0adv     = 1.0;          % J
P.bs        = [50 50];      % station de base

% Modele radio du premier ordre
P.k         = 4000;         % taille du paquet, bits
P.Eelec     = 50e-9;        % J/bit
P.efs       = 10e-12;       % J/bit/m^2
P.emp       = 0.0013e-12;   % J/bit/m^4
P.EDA       = 5e-9;         % J/bit/signal
P.d0        = sqrt(P.efs / P.emp);   % environ 87,7 m

% LEACH
P.p         = 0.05;         % probabilite d'etre tete de grappe

% PSO
P.K         = round(P.p * P.N);      % nombre de tetes de grappe par tour
P.swarm     = 30;           % particules
P.iters     = 100;          % iterations maximales
P.w         = 0.729;        % inertie
P.c1        = 1.494;
P.c2        = 1.494;
P.vmaxFrac  = 0.2;          % vitesse maximale, fraction de l'espace de recherche
P.weights   = [0.5 0.3 0.2];% alpha, beta, gamma
P.recluster = 5;            % reclustering tous les 5 tours

% Simulation
P.maxRounds = 5000;         % arret si des noeuds sont encore vivants
P.nRuns     = 50;           % executions independantes
end

function set_seed(s)
if exist('OCTAVE_VERSION', 'builtin')
    rand('twister', s);
    randn('state', s);
else
    rng(s);
end
end

function net = init_network(P)
% INIT_NETWORK  Deploiement aleatoire uniforme et energies initiales.
net.pos = P.area * rand(P.N, 2);
net.E0  = P.E0normal * ones(P.N, 1);
adv = randperm(P.N, P.nAdv);
net.E0(adv) = P.E0adv;
net.dBS = sqrt(sum((net.pos - repmat(P.bs, P.N, 1)).^2, 2));
end

function [Etx, Erx] = radio_energy(P, d)
% RADIO_ENERGY  Energie d'emission et de reception d'un paquet de P.k bits
% sur la distance d (vecteur), modele radio du premier ordre.
Etx = zeros(size(d));
fs = d < P.d0;
Etx(fs)  = P.k * P.Eelec + P.k * P.efs * d(fs).^2;
Etx(~fs) = P.k * P.Eelec + P.k * P.emp * d(~fs).^4;
Erx = P.k * P.Eelec;
end

function [cost, Davg, nPkt] = round_cost(P, net, alive, chIds)
% ROUND_COST  Energie consommee par chaque noeud pendant un tour.
%   Les membres envoient un paquet a leur tete la plus proche. Chaque tete
%   recoit, agrege m+1 signaux, puis envoie un paquet a la station de base.
%   Sans tete, chaque noeud vivant envoie directement a la station de base.
N = P.N;
cost = zeros(N, 1);
if isempty(chIds)
    d = net.dBS(alive);
    cost(alive) = radio_energy(P, d);
    Davg = mean(d);
    nPkt = nnz(alive);
    return
end
isCH = false(N, 1);
isCH(chIds) = true;
mem = find(alive & ~isCH);
K = numel(chIds);
if isempty(mem)
    cost(chIds) = radio_energy(P, net.dBS(chIds)) + P.k * P.EDA;
    Davg = 0;
    nPkt = K;
    return
end
dx = repmat(net.pos(mem, 1), 1, K) - repmat(net.pos(chIds, 1)', numel(mem), 1);
dy = repmat(net.pos(mem, 2), 1, K) - repmat(net.pos(chIds, 2)', numel(mem), 1);
D  = sqrt(dx.^2 + dy.^2);
[dmin, idx] = min(D, [], 2);
cost(mem) = radio_energy(P, dmin);
cnt = accumarray(idx, 1, [K 1]);
[~, Erx] = radio_energy(P, 0);
cost(chIds) = cnt * Erx + (cnt + 1) * P.k * P.EDA ...
              + radio_energy(P, net.dBS(chIds));
Davg = mean(dmin);
nPkt = K;
end

function [chIds, lastCH] = leach_select(r, P, alive, lastCH)
% LEACH_SELECT  Election LEACH standard.
%   Seuil T(n) = p / (1 - p * mod(r-1, 1/p)) pour les noeuds de l'ensemble G,
%   ceux qui n'ont pas encore ete tete durant l'epoque en cours. Une epoque
%   dure 1/p tours et G est reinitialise a chaque debut d'epoque. Les noeuds
%   morts ne sont jamais eligibles.
period = round(1 / P.p);
phase = mod(r - 1, period);
Tn = P.p / (1 - P.p * phase);
epochStart = r - phase;
eligible = alive & (lastCH < epochStart);
chosen = eligible & (rand(P.N, 1) < Tn);
chIds = find(chosen);
lastCH(chosen) = r;
end

function f = pso_fitness(P, net, E, alive, chIds)
% PSO_FITNESS  Fitness d'un jeu de tetes de grappe, a minimiser.
%   f = alpha * Davg + beta * (1 / Eres(CH)) + gamma * sigma_E
%   Chaque terme est normalise pour que les poids alpha, beta, gamma aient
%   un sens :
%     Davg    divise par la diagonale du champ
%     1/Eres  multiplie par l'energie initiale moyenne
%     sigma_E ecart type de l'energie residuelle PREDITE apres le tour,
%             divise par sa moyenne. L'ecart type de l'energie actuelle ne
%             depend pas du choix des tetes, il serait constant.
[cost, Davg] = round_cost(P, net, alive, chIds);
Epred = max(E - cost, 0);
Ea = Epred(alive);
Eref = mean(net.E0);
dRef = sqrt(2) * P.area;
t1 = Davg / dRef;
t2 = Eref / max(mean(E(chIds)), 1e-12);
t3 = std(Ea) / max(mean(Ea), 1e-12);
f = P.weights(1) * t1 + P.weights(2) * t2 + P.weights(3) * t3;
end

function chIds = pso_select(P, net, E, alive)
% PSO_SELECT  Choix des tetes de grappe par PSO.
%   Une particule est un vecteur continu de K positions dans [1, M], ou M
%   est le nombre de noeuds vivants. Elle est decodee en K indices distincts
%   de noeuds vivants, les doublons sont remplaces au hasard.
idxAlive = find(alive);
M = numel(idxAlive);
K = min(P.K, M);
if M <= K
    chIds = idxAlive;
    return
end
S = P.swarm;
x = 1 + (M - 1) * rand(S, K);
vmax = P.vmaxFrac * (M - 1);
v = (2 * rand(S, K) - 1) * vmax;
pBest = x;
pFit = inf(S, 1);
gFit = inf;
gBest = x(1, :);
for it = 1:P.iters
    for s = 1:S
        ids = decode_particle(x(s, :), M, K);
        f = pso_fitness(P, net, E, alive, idxAlive(ids));
        if f < pFit(s)
            pFit(s) = f;
            pBest(s, :) = x(s, :);
        end
        if f < gFit
            gFit = f;
            gBest = x(s, :);
        end
    end
    r1 = rand(S, K);
    r2 = rand(S, K);
    v = P.w * v + P.c1 * r1 .* (pBest - x) + P.c2 * r2 .* (repmat(gBest, S, 1) - x);
    v = max(min(v, vmax), -vmax);
    x = x + v;
    x = max(min(x, M), 1);
end
chIds = idxAlive(decode_particle(gBest, M, K));
end

function ids = decode_particle(x, M, K)
r = min(max(round(x), 1), M);
[~, first] = unique(r, 'first');
ids = r(sort(first));
if numel(ids) < K
    rest = setdiff(1:M, ids);
    ids = [ids, rest(randperm(numel(rest), K - numel(ids)))];
end
end

function res = run_protocol(proto, P, net)
% RUN_PROTOCOL  Simule 'leach' ou 'pso' sur le reseau net jusqu'a la mort de
% tous les noeuds, ou P.maxRounds. MAC ideal : pas de collision.
N = P.N;
E = net.E0;
alive = E > 0;
lastCH = -inf(N, 1);
chIds = [];
R = P.maxRounds;
res.aliveCnt = zeros(R, 1);
res.Eavg     = zeros(R, 1);
res.Evar     = zeros(R, 1);
res.Econs    = zeros(R, 1);
res.packets  = zeros(R, 1);
res.nCH      = zeros(R, 1);
res.FND = NaN;
res.LND = NaN;
res.rounds = 0;
for r = 1:R
    if ~any(alive)
        break
    end
    switch proto
        case 'leach'
            [chIds, lastCH] = leach_select(r, P, alive, lastCH);
        case 'pso'
            if mod(r - 1, P.recluster) == 0 || any(~alive(chIds))
                chIds = pso_select(P, net, E, alive);
            end
        otherwise
            error('protocole inconnu : %s', proto);
    end
    chIds = chIds(alive(chIds));
    [cost, ~, ~] = round_cost(P, net, alive, chIds);
    ok = E >= cost;
    delivered = ok(chIds);
    if isempty(chIds)
        nPkt = nnz(alive & ok);
    else
        nPkt = nnz(delivered);
    end
    spent = min(E, cost);
    E = E - spent;
    E(~ok) = 0;
    alive = E > 0;
    res.rounds = r;
    res.aliveCnt(r) = nnz(alive);
    res.Eavg(r)  = mean(E);
    res.Evar(r)  = var(E);
    res.Econs(r) = sum(spent);
    res.packets(r) = nPkt;
    res.nCH(r) = numel(chIds);
    if isnan(res.FND) && any(~alive)
        res.FND = r;
    end
    if ~any(alive)
        res.LND = r;
    end
end
end
