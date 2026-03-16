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
                       max(double(blk.rf_time_s)) * 1e3;
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

    % --- Plot (3 rows x 2 cols: col1=freq mod, col2=phase) -------------
    n_defs     = numel(ids);
    cmap       = lines(n_defs);
    bot_margin = 0.09;

    fig    = figure('Name', sprintf('Freqmod defs: %s', truth.base_name), 'Color', 'w');
    labels = {'Gx', 'Gy', 'Gz'};
    fm_ax  = gobjects(1, 3);
    ph_ax  = gobjects(1, 3);

    for ax = 1:3
        fm_ax(ax) = subplot(3, 2, 2*(ax-1) + 1);
        hold(fm_ax(ax), 'on');
        ylabel(fm_ax(ax), sprintf('%s (Hz/m)', labels{ax}));
        grid(fm_ax(ax), 'on');
        if ax == 1, title(fm_ax(ax), 'Freq Mod'); end

        ph_ax(ax) = subplot(3, 2, 2*(ax-1) + 2);
        hold(ph_ax(ax), 'on');
        ylabel(ph_ax(ax), sprintf('\\phi_{%s} (rad/m)', labels{ax}));
        grid(ph_ax(ax), 'on');
        if ax == 1, title(ph_ax(ax), 'Phase'); end
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

        t0_ms    = def_t0_ms(i);
        dur_ms   = double(def.duration_us) * 1e-3;
        t_win_ms = (0:double(def.num_samples)-1).' * double(def.raster_us) * 1e-3;

        for ax = 1:3
            fm = double(def.waveform(:, ax));
            ph = cumsum(fm * 2 * pi * dt_s);

            % Zero-padded envelope: 0 -> step up at t0 -> waveform ->
            % step down at t0+duration -> 0 until end of segment.
            t_plot  = [0;  t0_ms;  t0_ms + t_win_ms;  t0_ms + dur_ms;  total_dur_ms];
            fm_plot = [0;  0;      fm;                 0;               0            ];
            ph_plot = [0;  0;      ph;                 ph(end);         ph(end)      ];

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

    sgtitle(fig, sprintf('Frequency Modulation Definitions (%s)  seg %d', ...
        truth.base_name, s_idx - 1), 'Interpreter', 'none');
end

function truth = coerce_truth(base_or_truth)
    if isstruct(base_or_truth)
        truth = base_or_truth;
    else
        truth = testutils.truth_parse_case(base_or_truth);
    end
end
