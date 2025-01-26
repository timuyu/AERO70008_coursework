%% Admin
close all;
clc;
clear;

%% Question 1

% scheme = "LeapFrog";
scheme = "AngledDerivative";

% Set-up
xd = [-1, 1];
Td = [0, 1];

a = 0.75;
c = 0.6;

dxs = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001];
dts = c * dxs / a;

% Iterate
fig = figure(1);
    labs = "x = " + dxs;

for ix = 1:length(dxs)
    dx = dxs(ix);
    xs = xd(1):dx:xd(2);
    nx = length(xs);

    dt = dts(ix);
    ts = Td(1):dt:Td(2);
    nt = length(ts);

    indt = [1, find(abs(ts - 0.5) == min(abs(ts - 0.5)), 1), find(abs(ts - 1.0) == min(abs(ts - 1.0)), 1)];

    u = Solver1(scheme, nt, nx, ts, xs, a, c);

    Save1(scheme, ix, dx, dt, a, c, ts, xs, u)

    for i = 1:length(indt)
        fig = Plot1(fig, i, xs, u, indt, ts, xd, dx, dxs);
    end

end
hold off;

% Error
eps = zeros(length(dxs), 2);

for ischeme = 1:2
    for igrid = 1:length(dxs)

        if ischeme == 1, fin = sprintf("CW1_data/Q1/LF_dx%i", igrid); end
        if ischeme == 2, fin = sprintf("CW1_data/Q1/AD_dx%i", igrid); end

        load(fin, 'xs', 'u')
        ut = u(1:end-1, end);

        xs = xs(1:end-1);

        u_exact = zeros(length(xs), 1);
        u_exact = Set_IC1(u_exact, xs, 1, a);

        % eps(igrid, ischeme) = max(abs(u_exact - ut));
        eps(igrid, ischeme) = norm(u_exact - ut, 2);

    end
end

figure(2)
    loglog(dxs, eps(:, 1), 'Color', [0      0.4470 0.7410], 'LineWidth', 1); hold on;
    loglog(dxs, eps(:, 2), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1); hold off;
        xlabel('$\Delta x$', 'Interpreter', 'latex')
        ylabel('$\epsilon$', 'Interpreter', 'latex')
        legend(["Leap Frog", "Angled Derivative"], 'Interpreter', 'latex')
        grid on
        set(gca, 'FontSize', 18)


%% Question 2
clear;

schemes = ["LeapFrog", "AngledDerivative"];
scheme_titles = ["Leap Frog", "Angled Derivative"];

xd = [-1, 1];
Td = [0, 0.5];

a = 0.75;
c = 0.6;
dx = 0.005;

dt = c * dx / a;

xs = xd(1):dx:xd(2);
ts = Td(1):dt:Td(2);

nx = length(xs);
nt = length(ts);

fig = figure(4);

for ischeme = 1:length(schemes)
    scheme = schemes(ischeme);
    scheme_title = scheme_titles(ischeme);

    u = Solver2(scheme, nt, nx, xs, c);
    
    Save2(scheme, dx, dt, a, c, ts, xs, u)
    
    if ischeme == 1, fin = sprintf("CW1_data/Q2/LF"); end
    if ischeme == 2, fin = sprintf("CW1_data/Q2/AD"); end

    load(fin, 'xs', 'u')
    ut = u(:, end);

    fig = Plot2(fig, ischeme, scheme_title, xs, ts, u, a, nt, xd);
end


%% Functions
function [u] = Set_IC1(u, x, t, a)  % I recognise that MATLAB doesn't Pass by Reference so no need to pass u here, but I like to do it to be reminded

    u(:) = sin(2*pi*(x - a*t));

end

function [u] = Set_IC2(u, x, i, c, n)

    if i == 1
        u((x >= -0.2) & (x <= 0.2)) = 1;
    else
        u = LaxWendroff(n, u, c);
    end

end

function [u3] = LeapFrog(n, u1, u2, u3, c)  % I recognise that u3 isn't Passed by
                                            % Reference, but I like to do it to be reminded

    u3(1)     = u1(1)     - c * ( u2(2)     - u2(n-1));
    u3(2:n-2) = u1(2:n-2) - c * ( u2(3:n-1) - u2(1:n-3) );                  % interior
    u3(n-1)   = u1(n-1)   - c * ( u2(1)     - u2(n-2));

    u3(n) = u3(1);                                                          % periodic bcs

end

function [u3] = AngledDerivative(n, u1, u2, u3, c)

    u3(1)     = u1(n-1)   + (1 - 2*c)*( u2(1)     - u2(n-1));
    u3(2:n-1) = u1(1:n-2) + (1 - 2*c)*( u2(2:n-1) - u2(1:n-2));             % interior

    u3(n) = u3(1);                                                          % periodic bcs

end

function [u2] = LaxWendroff(n, u1, c)

    u2(1)     = 0.5 * c*(1 + c) * u1(n-1)   + (1 - c^2) * u1(1)     - 0.5*c*(1 - c) * u1(2);
    u2(2:n-2) = 0.5 * c*(1 + c) * u1(1:n-3) + (1 - c^2) * u1(2:n-2) - 0.5*c*(1 - c) * u1(3:n-1);
    u2(n-1)   = 0.5 * c*(1 + c) * u1(n-2)   + (1 - c^2) * u1(n-1)   - 0.5*c*(1 - c) * u1(1);

    u2(n) = u2(1);

end

function [u] = Solver1(scheme, nt, nx, ts, xs, a, c)

    u = zeros(nx, nt);
    ut = zeros(nx, 1);

    for i = 1:nt
        if (i == 1) || (i == 2)
            ut = Set_IC1(ut, xs, ts(i), a);     % sets ics using exact solution

        else
            if scheme == "LeapFrog"
                ut = LeapFrog(nx, u(:, i-2), u(:, i-1), ut, c);

            elseif scheme == "AngledDerivative"
                ut = AngledDerivative(nx, u(:, i-2), u(:, i-1), ut, c);

            end
        end
        
        u(:, i) = ut;
    end
end

function [] = Save1(scheme, ix, dx, dt, a, c, ts, xs, u)

    if scheme == "LeapFrog"
        save(sprintf("CW1_data/Q1/LF_dx%i", ix), 'dx', 'dt', 'a', 'c', 'ts', 'xs', 'u')
    elseif scheme == "AngledDerivative"
        save(sprintf("CW1_data/Q1/AD_dx%i", ix), 'dx', 'dt', 'a', 'c', 'ts', 'xs', 'u')
    end

end

function [fig] = Plot1(fig, i, xs, u, indt, ts, xd, dx, dxs)

    figure(fig)
        % subplot(3, 2, i)
            subplot(3, 1, i)
            plot(xs, u(:, indt(i)), 'LineWidth', 1); hold on;
                xlim(xd);
                ylim(xd);
        if dx == dxs(end)
            xlabel('$x$', 'Interpreter', 'latex')
            ylabel('$u$', 'Interpreter', 'latex')
            title(sprintf('$t = %1.2f$', ts(indt(i))), 'Interpreter', 'latex')
            set(gca, 'FontSize', 18)
        end
        % if i == length(indt) && ix == length(dxs)
        %     ax = subplot(3, 2, 6, 'Visible', 'off');
        %     ax = subplot(3, 2, 6, 'Visible', 'off');
        %     axPos = ax.Position; delete(ax)
        % 
        %     hl = legend(["$\Delta x = 0.1$", "$\Delta x = 0.05$", "$\Delta x = 0.02$", "$\Delta x = 0.01$", "$\Delta x = 0.005$", "$\Delta x = 0.002$", "$\Delta x = 0.00$1"], 'Interpreter', 'latex', 'FontSize', 12);
        %     % set(gca, 'FontSize', 12)
        %     hl.Position(1:2) = axPos(1:2);
        % end
        grid on

end

function [u] = Solver2(scheme, nt, nx, xs, c)

    u = zeros(nx, nt);
    ut = zeros(nx, 1);

    for i = 1:nt
        
        if (i == 1) || (i == 2)
            ut = Set_IC2(ut, xs, i, c, nx);                         % sets ics using exact solution

        else
            if scheme == "LeapFrog"
                ut = LeapFrog(nx, u(:, i-2), u(:, i-1), ut, c);        % Uses Leap Frog numerical scheme

            elseif scheme == "AngledDerivative"
                ut = AngledDerivative(nx, u(:, i-2), u(:, i-1), ut, c);

            end
        end

        u(:, i) = ut;
    end
end

function [] = Save2(scheme, dx, dt, a, c, ts, xs, u)

    if scheme == "LeapFrog"
        save(sprintf("CW1_data/Q2/LF"), 'dx', 'dt', 'a', 'c', 'ts', 'xs', 'u')
    elseif scheme == "AngledDerivative"
        save(sprintf("CW1_data/Q2/AD"), 'dx', 'dt', 'a', 'c', 'ts', 'xs', 'u')
    end

end

function [fig] = Plot2(fig, ischeme, scheme, xs, ts, u, a, it, xd)

    figure(fig)
    subplot(2, 1, ischeme)
        plot(xs + a*ts(it), u(:, 1), '--k', 'LineWidth', 0.5); hold on;
        plot(xs, u(:, it), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1); hold off;
            xlim(xd);
            ylim(xd + 0.5);
            xlabel('$x$', 'Interpreter', 'latex')
            ylabel('$u$', 'Interpreter', 'latex')
            title(sprintf("%s", scheme), 'Interpreter', 'latex')
            grid on;
            set(gca, 'FontSize', 18)
end