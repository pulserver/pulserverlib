function fig = truth_plot_freqmod_defs(base_or_truth, varargin)
%TRUTH_PLOT_FREQMOD_DEFS Plot frequency modulation waveforms and cumulative phase.
%
% Rows 1-3: freq mod waveform per axis in Hz/m.
% Rows 4-6: cumulative phase integral 2pi*integral(G dt) in rad/m.
% A circle marker on the phase rows indicates the stored ref_integral value
% at ref_time_us as a sanity check of the running integral.
% All six rows share the same x-axis in us.
%
%   truth_plot_freqmod_defs('gre_2d_1sl_1avg')
%   truth_plot_freqmod_defs(truth_struct, 'def_idx', [1 2])

    p = inputParser;
    addParameter(p, 'def_idx', [], @(x) isempty(x) || isnumeric(x));
    parse(p, varargin{:});

    truth = coerce_truth(base_or_truth);

    if truth.freqmod_def.num_defs == 0
        error('truth:empty', 'No frequency modulation definitions found for %s', truth.base_name);
    end

    if isempty(p.Results.def_idx)
        ids = 1:truth.freqmod_def.num_defs;
    else
        ids = p.Results.def_idx(:).';
    end

    n_defs  = numel(ids);
    cmap    = lines(n_defs);
    bot_margin = 0.09;

    fig = figure('Name', sprintf('Freqmod defs: %s', truth.base_name), 'Color', 'w');

    labels = {'Gx', 'Gy', 'Gz'};
    fm_ax  = gobjects(1, 3);
    ph_ax  = gobjects(1, 3);

    for ax = 1:3
        fm_ax(ax) = subplot(6, 1, ax);
        hold(fm_ax(ax), 'on');
        ylabel(fm_ax(ax), sprintf('%s (Hz/m)', labels{ax}));
        grid(fm_ax(ax), 'on');

        ph_ax(ax) = subplot(6, 1, 3 + ax);
        hold(ph_ax(ax), 'on');
        ylabel(ph_ax(ax), sprintf('\\phi_{%s} (rad/m)', labels{ax}));
        grid(ph_ax(ax), 'on');
    end
    xlabel(fm_ax(3), 'Time (us)');
    xlabel(ph_ax(3), 'Time (us)');

    def_lbls = cell(1, n_defs);
    for i = 1:n_defs
        d = ids(i);
        if d < 1 || d > truth.freqmod_def.num_defs
            error('truth:input', 'def_idx %d out of range', d);
        end
        def    = truth.freqmod_def.defs(d);
        t_us   = (0:double(def.num_samples)-1).' * double(def.raster_us);
        dt_s   = double(def.raster_us) * 1e-6;
        kind   = 'RF';
        if def.type == 1, kind = 'ADC'; end
        def_lbls{i} = sprintf('Def %d (%s)', d, kind);

        for ax = 1:3
            fm  = double(def.waveform(:, ax));
            ph  = cumsum(fm * 2 * pi * dt_s);   % rad/m

            plot(fm_ax(ax), t_us, fm, 'Color', cmap(i,:), 'LineWidth', 1.2);
            plot(ph_ax(ax), t_us, ph, 'Color', cmap(i,:), 'LineWidth', 1.2);

            % Stored ref_integral (Hz*s/m) * 2pi = rad/m; marker at ref_time_us.
            ref_ph = double(def.ref_integral(ax)) * 2 * pi;
            plot(ph_ax(ax), double(def.ref_time_us), ref_ph, 'o', ...
                'Color', cmap(i,:), 'MarkerSize', 6, 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
    end

    all_ax = [fm_ax, ph_ax];
    linkaxes(all_ax, 'x');

    % Compress subplots upward to leave room for the global legend.
    for k = 1:6
        pos    = get(all_ax(k), 'Position');
        pos(2) = pos(2) * (1 - bot_margin) + bot_margin;
        pos(4) = pos(4) * (1 - bot_margin);
        set(all_ax(k), 'Position', pos);
    end

    % Global horizontal legend in a hidden axes at the bottom.
    lgd_ax = axes('Parent', fig, ...
        'Position', [0.08, 0.005, 0.84, bot_margin - 0.01], ...
        'Visible', 'off');
    hold(lgd_ax, 'on');
    dummy_h = gobjects(1, n_defs);
    for i = 1:n_defs
        dummy_h(i) = plot(lgd_ax, NaN, NaN, '-', 'Color', cmap(i,:), 'LineWidth', 1.2);
    end
    legend(lgd_ax, dummy_h, def_lbls, ...
        'Orientation', 'horizontal', 'Location', 'north', 'Box', 'off');

    sgtitle(sprintf('Frequency Modulation Definitions (%s)', truth.base_name), 'Interpreter', 'none');
end

function truth = coerce_truth(base_or_truth)
    if isstruct(base_or_truth)
        truth = base_or_truth;
    else
        truth = testutils.truth_parse_case(base_or_truth);
    end
end
%
% Rows 1-3 show the freq mod waveform per axis in Hz/m.
% Rows 4-6 show the cumulative phase integral (2pi * int G dt) in rad/m.
% A circle marker on the phase plots indicates the stored ref_integral
% value at ref_time_us, providing a quick sanity-check of the integral.
%
%   truth_plot_freqmod_defs('gre_2d_1sl_1avg')
%   truth_plot_freqmod_defs(truth_struct, 'def_idx', [1 2])

    p = inputParser;
    addParameter(p, 'def_idx', [], @(x) isempty(x) || isnumeric(x));
    parse(p, varargin{:});

    truth = coerce_truth(base_or_truth);

    if truth.freqmod_def.num_defs == 0
        error('truth:empty', 'No frequency modulation definitions found for %s', truth.base_name);
    end

    if isempty(p.Results.def_idx)
        ids = 1:truth.freqmod_def.num_defs;
    else
        ids = p.Results.def_idx(:).';
    end

    fig = figure('Name', sprintf('Freqmod defs: %s', truth.base_name), 'Color', 'w');

    labels = {'Gx', 'Gy', 'Gz'};
    fm_ax = gobjects(1, 3);   % freq mod (Hz/m)
    ph_ax = gobjects(1, 3);   % cumulative phase (rad/m)

    for ax = 1:3
        fm_ax(ax) = subplot(6, 1, ax);
        hold(fm_ax(ax), 'on');
        ylabel(fm_ax(ax), sprintf('%s (Hz/m)', labels{ax}));
        grid(fm_ax(ax), 'on');

        ph_ax(ax) = subplot(6, 1, 3 + ax);
        hold(ph_ax(ax), 'on');
        ylabel(ph_ax(ax), sprintf('\\phi_{%s} (rad/m)', labels{ax}));
        grid(ph_ax(ax), 'on');
    end
    xlabel(fm_ax(3), 'Time (us)');
    xlabel(ph_ax(3), 'Time (us)');

    lgd = {};
    for i = 1:numel(ids)
        d = ids(i);
        if d < 1 || d > truth.freqmod_def.num_defs
            error('truth:input', 'def_idx %d out of range', d);
        end
        def = truth.freqmod_def.defs(d);
        t = (0:double(def.num_samples)-1).' * double(def.raster_us);
        dt_s = double(def.raster_us) * 1e-6;

        kind = 'RF';
        if def.type == 1
            kind = 'ADC';
        end
        lgd{end + 1} = sprintf('Def %d (%s)', d, kind); %#ok<AGROW>

        for ax = 1:3
            fm = double(def.waveform(:, ax));
            ph = cumsum(fm * 2 * pi * dt_s);   % rad/m

            plot(fm_ax(ax), t, fm, 'LineWidth', 1.2);
            plot(ph_ax(ax), t, ph, 'LineWidth', 1.2);

            % Mark the stored reference integral at ref_time_us.
            % ref_integral is in Hz*s/m; multiply by 2pi to get rad/m.
            ref_ph = double(def.ref_integral(ax)) * 2 * pi;
            plot(ph_ax(ax), double(def.ref_time_us), ref_ph, 'o', ...
                'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end

    legend(fm_ax(1), lgd, 'Location', 'eastoutside');
    linkaxes([fm_ax, ph_ax], 'x');
    sgtitle(sprintf('Frequency Modulation Definitions (%s)', truth.base_name), 'Interpreter', 'none');
end

function truth = coerce_truth(base_or_truth)
    if isstruct(base_or_truth)
        truth = base_or_truth;
    else
        truth = testutils.truth_parse_case(base_or_truth);
    end
end
