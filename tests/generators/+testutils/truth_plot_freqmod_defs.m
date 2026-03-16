function fig = truth_plot_freqmod_defs(base_or_truth, varargin)
%TRUTH_PLOT_FREQMOD_DEFS Plot frequency modulation waveforms and cumulative phase.
%
% Each def is shown on the full segment timeline in ms: zero before the
% active window (RF or ADC delay), the stored waveform during the window,
% and zero again after.  Rows 1-3 show Hz/m; rows 4-6 show rad/m.
% A circle marker on the phase rows checks the stored ref_integral value.
%
%   truth_plot_freqmod_defs('gre_2d_1sl_1avg')
%   truth_plot_freqmod_defs(truth_struct, 'def_idx', [1 2])
%   truth_plot_freqmod_defs(truth_struct, 'segment_idx', 2)

    p = inputParser;
    addParameter(p, 'def_idx',     [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'segment_idx', 1,  @(x) isnumeric(x) && isscalar(x) && x >= 1);
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

    adc_dwell_ms = double(truth.meta.adc_dwell_ns) * 1e-6;
    n_adc        = double(truth.meta.adc_samples);
    adc_dur_ms   = n_adc * adc_dwell_ms;

    % --- Build segment timing context (ms) ----------------------------
    s_idx = p.Results.segment_idx;
    if s_idx < 1 || s_idx > truth.segment_def.num_segments
        error('truth:input', 'segment_idx %d out of range', s_idx);
    end
    seg      = truth.segment_def.segments(s_idx);
    n_blocks = seg.num_blocks;

    block_end_ms = zeros(1, n_blocks);
    for b = 1:n_blocks
        blk   = seg.blocks(b);
        spans = zeros(1, 5);
        for ch = 1:3
            n = double(blk.grad_n(ch));
            if n > 0
                spans(ch) = double(blk.grad_delay(ch)) * 1e3 + ...
                            max(double(blk.grad_time_s{ch})) * 1e3;
            end
        end
        rf_n = double(blk.rf_n);
        if rf_n > 0
            spans(4) = double(blk.rf_delay) * 1e3 + ...
                       (rf_n - 1) * double(blk.rf_raster_us) * 1e-3;
        end
        if blk.has_adc
            spans(5) = double(blk.adc_delay) * 1e3 + adc_dur_ms;
        end
        block_end_ms(b) = max(spans);
    end
    block_start_ms = [0, cumsum(block_end_ms + 0.002)];
    total_dur_ms   = block_start_ms(end);

    % --- Match each def to its source block (RF or ADC type) ----------
    rf_fm_blks  = [];
    adc_fm_blks = [];
    for b = 1:n_blocks
        blk = seg.blocks(b);
        if ~blk.has_freq_mod, continue; end
        if blk.has_rf,  rf_fm_blks(end+1)  = b; end %#ok<AGROW>
        if blk.has_adc, adc_fm_blks(end+1) = b; end %#ok<AGROW>
    end

    rf_ctr  = 0;
    adc_ctr = 0;
    def_t0_ms = zeros(1, numel(ids));   % active-window start (ms)
    for i = 1:numel(ids)
        def = truth.freqmod_def.defs(ids(i));
        if def.type == 0
            rf_ctr = rf_ctr + 1;
            if rf_ctr <= numel(rf_fm_blks)
                b = rf_fm_blks(rf_ctr);
                def_t0_ms(i) = block_start_ms(b) + double(seg.blocks(b).rf_delay) * 1e3;
            end
        else
            adc_ctr = adc_ctr + 1;
            if adc_ctr <= numel(adc_fm_blks)
                b = adc_fm_blks(adc_ctr);
                def_t0_ms(i) = block_start_ms(b) + double(seg.blocks(b).adc_delay) * 1e3;
            end
        end
    end

    % --- Plot ---------------------------------------------------------
    n_defs     = numel(ids);
    cmap       = lines(n_defs);
    bot_margin = 0.09;

    fig    = figure('Name', sprintf('Freqmod defs: %s', truth.base_name), 'Color', 'w');
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
    xlabel(fm_ax(3), 'Time (ms)');
    xlabel(ph_ax(3), 'Time (ms)');

    def_lbls = cell(1, n_defs);
    for i = 1:n_defs
        d    = ids(i);
        def  = truth.freqmod_def.defs(d);
        dt_s = double(def.raster_us) * 1e-6;
        kind = 'RF';
        if def.type == 1, kind = 'ADC'; end
        def_lbls{i} = sprintf('Def %d (%s)', d, kind);

        t0_ms   = def_t0_ms(i);
        t_win_ms = (0:double(def.num_samples)-1).' * double(def.raster_us) * 1e-3;

        for ax = 1:3
            fm = double(def.waveform(:, ax));
            ph = cumsum(fm * 2 * pi * dt_s);

            % Zero-padded envelope spanning the full segment timeline.
            t_plot  = [0;    t0_ms;  t0_ms + t_win_ms;  total_dur_ms];
            fm_plot = [0;    0;      fm;                 0           ];
            ph_plot = [0;    0;      ph;                 ph(end)     ];

            plot(fm_ax(ax), t_plot, fm_plot, 'Color', cmap(i,:), 'LineWidth', 1.2);
            plot(ph_ax(ax), t_plot, ph_plot, 'Color', cmap(i,:), 'LineWidth', 1.2);

            % ref_integral (Hz·s/m) × 2π = rad/m; marker at ref_time.
            ref_ph = double(def.ref_integral(ax)) * 2 * pi;
            plot(ph_ax(ax), t0_ms + double(def.ref_time_us) * 1e-3, ref_ph, 'o', ...
                'Color', cmap(i,:), 'MarkerSize', 6, 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
    end

    all_ax = [fm_ax, ph_ax];
    linkaxes(all_ax, 'x');

    for k = 1:6
        pos    = get(all_ax(k), 'Position');
        pos(2) = pos(2) * (1 - bot_margin) + bot_margin;
        pos(4) = pos(4) * (1 - bot_margin);
        set(all_ax(k), 'Position', pos);
    end

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

    sgtitle(sprintf('Frequency Modulation Definitions (%s)  seg %d', ...
        truth.base_name, s_idx), 'Interpreter', 'none');
end
%
% Each def is shown on the full segment timeline: zero before the active
% window (RF or ADC delay), the stored waveform during the window, and
% zero again after.  Rows 1-3 show Hz/m; rows 4-6 show rad/m (cumulative
% phase integral 2*pi*integral(G dt), held flat after the window).
% A circle marker on the phase rows indicates the stored ref_integral
% value at ref_time_us as a sanity check of the running integral.
%
% The segment used for timing context is selected by 'segment_idx' (1).
%   def type==0 (RF)  matches the N-th block with has_rf+has_freq_mod.
%   def type==1 (ADC) matches the N-th block with has_adc+has_freq_mod.
%
%   truth_plot_freqmod_defs('gre_2d_1sl_1avg')
%   truth_plot_freqmod_defs(truth_struct, 'def_idx', [1 2])
%   truth_plot_freqmod_defs(truth_struct, 'segment_idx', 2)

    p = inputParser;
    addParameter(p, 'def_idx',       [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'segment_idx',   1,  @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'grad_raster_us', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
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

    % --- Derive gradient raster ----------------------------------------
    grad_raster_us = p.Results.grad_raster_us;
    if grad_raster_us <= 0
        grad_raster_us = double(truth.freqmod_def.defs(1).raster_us);
    end
    adc_dwell_us = double(truth.meta.adc_dwell_ns) / 1000;
    n_adc        = double(truth.meta.adc_samples);
    adc_dur_us   = n_adc * adc_dwell_us;

    % --- Build segment timing context ----------------------------------
    s_idx = p.Results.segment_idx;
    if s_idx < 1 || s_idx > truth.segment_def.num_segments
        error('truth:input', 'segment_idx %d out of range', s_idx);
    end
    seg      = truth.segment_def.segments(s_idx);
    n_blocks = seg.num_blocks;

    block_width_us = zeros(1, n_blocks);
    for b = 1:n_blocks
        blk   = seg.blocks(b);
        spans = zeros(1, 5);
        for ch = 1:3
            spans(ch) = double(blk.grad_delay(ch)) * 1e6 + double(blk.grad_n(ch)) * grad_raster_us;
        end
        spans(4) = double(blk.rf_delay) * 1e6 + double(blk.rf_n) * grad_raster_us;
        if blk.has_adc
            spans(5) = double(blk.adc_delay) * 1e6 + adc_dur_us;
        end
        block_width_us(b) = max(spans);
    end
    block_start_us = [0, cumsum(block_width_us + 2)];
    total_dur_us   = block_start_us(end);

    % --- Match each def to segment blocks (RF or ADC type) -------------
    rf_fm_blks  = [];
    adc_fm_blks = [];
    for b = 1:n_blocks
        blk = seg.blocks(b);
        if ~blk.has_freq_mod, continue; end
        if blk.has_rf,  rf_fm_blks(end+1)  = b; end %#ok<AGROW>
        if blk.has_adc, adc_fm_blks(end+1) = b; end %#ok<AGROW>
    end

    rf_ctr  = 0;
    adc_ctr = 0;
    def_t0  = zeros(1, numel(ids));   % active-window start (us) in segment timeline
    for i = 1:numel(ids)
        def = truth.freqmod_def.defs(ids(i));
        if def.type == 0
            rf_ctr = rf_ctr + 1;
            if rf_ctr <= numel(rf_fm_blks)
                b = rf_fm_blks(rf_ctr);
                def_t0(i) = block_start_us(b) + double(seg.blocks(b).rf_delay) * 1e6;
            end
        else
            adc_ctr = adc_ctr + 1;
            if adc_ctr <= numel(adc_fm_blks)
                b = adc_fm_blks(adc_ctr);
                def_t0(i) = block_start_us(b) + double(seg.blocks(b).adc_delay) * 1e6;
            end
        end
    end

    % --- Plot ----------------------------------------------------------
    n_defs     = numel(ids);
    cmap       = lines(n_defs);
    bot_margin = 0.09;

    fig    = figure('Name', sprintf('Freqmod defs: %s', truth.base_name), 'Color', 'w');
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
        d   = ids(i);
        def = truth.freqmod_def.defs(d);
        dt_s = double(def.raster_us) * 1e-6;
        kind = 'RF';
        if def.type == 1, kind = 'ADC'; end
        def_lbls{i} = sprintf('Def %d (%s)', d, kind);

        t0     = def_t0(i);
        t_win  = (0:double(def.num_samples)-1).' * double(def.raster_us);  % local time

        for ax = 1:3
            fm = double(def.waveform(:, ax));
            ph = cumsum(fm * 2 * pi * dt_s);

            % Zero-padded time + value vectors spanning [0 .. total_dur_us].
            % Step up at t0, active waveform, step down to 0 after window.
            t_plot  = [0;   t0;  t0 + t_win;          total_dur_us];
            fm_plot = [0;   0;   fm;                   0           ];
            ph_plot = [0;   0;   ph;                   ph(end)     ];

            plot(fm_ax(ax), t_plot, fm_plot, 'Color', cmap(i,:), 'LineWidth', 1.2);
            plot(ph_ax(ax), t_plot, ph_plot, 'Color', cmap(i,:), 'LineWidth', 1.2);

            % Ref-integral marker (ref_integral in Hz*s/m → *2pi = rad/m).
            ref_ph = double(def.ref_integral(ax)) * 2 * pi;
            plot(ph_ax(ax), t0 + double(def.ref_time_us), ref_ph, 'o', ...
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

    % Global horizontal legend at the bottom.
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

    sgtitle(sprintf('Frequency Modulation Definitions (%s)  seg %d', ...
        truth.base_name, s_idx), 'Interpreter', 'none');
end

function truth = coerce_truth(base_or_truth)
    if isstruct(base_or_truth)
        truth = base_or_truth;
    else
        truth = testutils.truth_parse_case(base_or_truth);
    end
end
