function figs = truth_plot_segments(base_or_truth, varargin)
%TRUTH_PLOT_SEGMENTS Plot segment block waveforms on a shared real-time axis.
%
% All five channels (Gx/Gy/Gz/RF/ADC) share the same x-axis in us.
% Block boundary positions are computed from the maximum event span per
% block given the gradient raster (auto-detected from freq-mod defs or TR
% waveforms). Gradients in mT/m, RF |B1| in uT, ADC as a binary mask.
%
% Note: block delays are stored in seconds (Pulseq convention) and are
% converted to us here. Waveform samples are scaled by grad_raster_us.
% Trap gradient keypoints (3-4 pts) are plotted at the same spacing as
% arbitrary waveform samples (approximate timing within the block).
%
%   truth_plot_segments('gre_2d_1sl_1avg')
%   truth_plot_segments(truth_struct, 'segment_idx', [1 2])
%   truth_plot_segments(truth_struct, 'grad_raster_us', 10)

    p = inputParser;
    addParameter(p, 'segment_idx', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'grad_raster_us', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, varargin{:});

    truth = coerce_truth(base_or_truth);

    GAMMA_HZ_PER_T = 42.577e6;
    GAMMA_GRAD     = GAMMA_HZ_PER_T * 1e-3;   % Hz/m -> mT/m: divide by this
    GAMMA_RF       = GAMMA_HZ_PER_T * 1e-6;   % Hz   -> uT:   divide by this

    % Derive gradient raster (samples -> us).
    grad_raster_us = p.Results.grad_raster_us;
    if grad_raster_us <= 0
        if truth.freqmod_def.num_defs > 0
            grad_raster_us = double(truth.freqmod_def.defs(1).raster_us);
        elseif truth.tr_waveforms.num_trs > 0 && numel(truth.tr_waveforms.waveforms(1).time_us) >= 2
            % TR waveforms are at 0.5 * gradRasterTime; segment waveforms at gradRasterTime.
            tr_step = double(truth.tr_waveforms.waveforms(1).time_us(2)) - ...
                      double(truth.tr_waveforms.waveforms(1).time_us(1));
            grad_raster_us = 2 * tr_step;
        else
            grad_raster_us = double(truth.meta.adc_dwell_ns) / 1000;
        end
    end

    adc_dwell_us = double(truth.meta.adc_dwell_ns) / 1000;
    n_adc        = double(truth.meta.adc_samples);
    adc_dur_us   = n_adc * adc_dwell_us;

    if isempty(p.Results.segment_idx)
        seg_ids = 1:truth.segment_def.num_segments;
    else
        seg_ids = p.Results.segment_idx(:).';
    end

    figs        = gobjects(1, numel(seg_ids));
    bot_margin  = 0.09;   % fraction of figure height reserved for the legend

    for ii = 1:numel(seg_ids)
        s = seg_ids(ii);
        if s < 1 || s > truth.segment_def.num_segments
            error('truth:input', 'segment_idx %d out of range', s);
        end

        seg      = truth.segment_def.segments(s);
        n_blocks = seg.num_blocks;
        cmap     = lines(n_blocks);

        fig  = figure('Name', sprintf('Segment %d: %s', s - 1, truth.base_name), 'Color', 'w');
        ax_h = gobjects(1, 5);

        % Block start positions in us (delays are stored in seconds -> *1e6).
        block_width_us = zeros(1, n_blocks);
        for b = 1:n_blocks
            blk    = seg.blocks(b);
            spans  = zeros(1, 5);
            for ch = 1:3
                d_us      = double(blk.grad_delay(ch)) * 1e6;
                spans(ch) = d_us + double(blk.grad_n(ch)) * grad_raster_us;
            end
            spans(4) = double(blk.rf_delay) * 1e6 + double(blk.rf_n) * grad_raster_us;
            if blk.has_adc
                spans(5) = double(blk.adc_delay) * 1e6 + adc_dur_us;
            end
            block_width_us(b) = max(spans);
        end
        block_start_us = [0, cumsum(block_width_us + 2)];  % 2 us visual gap between blocks

        % --- Gradient channels ---
        grad_ylabels = {'Gx (mT/m)', 'Gy (mT/m)', 'Gz (mT/m)'};
        for ch = 1:3
            ax_h(ch) = subplot(5, 1, ch);
            hold on;
            for b = 1:n_blocks
                blk = seg.blocks(b);
                n   = double(blk.grad_n(ch));
                if n > 0
                    y    = double(blk.grad_amp(ch)) * double(blk.grad_wave{ch}(:)) / GAMMA_GRAD;
                    t_us = block_start_us(b) + double(blk.grad_delay(ch)) * 1e6 + (0:n-1).' * grad_raster_us;
                    plot(t_us, y, 'Color', cmap(b,:), 'LineWidth', 1.1);
                end
            end
            hold off;
            ylabel(grad_ylabels{ch});
            grid on;
        end

        % --- RF |B1| ---
        ax_h(4) = subplot(5, 1, 4);
        hold on;
        for b = 1:n_blocks
            blk = seg.blocks(b);
            n   = double(blk.rf_n);
            if n > 0
                y    = double(blk.rf_amp) * double(blk.rf_rho(:)) / GAMMA_RF;
                t_us = block_start_us(b) + double(blk.rf_delay) * 1e6 + (0:n-1).' * grad_raster_us;
                plot(t_us, y, 'Color', cmap(b,:), 'LineWidth', 1.1);
            end
        end
        hold off;
        ylabel('RF |B1| (uT)');
        grid on;

        % --- ADC mask ---
        ax_h(5) = subplot(5, 1, 5);
        hold on;
        for b = 1:n_blocks
            blk = seg.blocks(b);
            if blk.has_adc && n_adc > 0
                t0_us = block_start_us(b) + double(blk.adc_delay) * 1e6;
                t1_us = t0_us + adc_dur_us;
                fill([t0_us, t1_us, t1_us, t0_us], [0, 0, 1, 1], cmap(b,:), ...
                    'FaceAlpha', 0.45, 'EdgeColor', 'none');
            end
        end
        hold off;
        ylim([-0.1, 1.5]);
        yticks([0, 1]);
        yticklabels({'0', '1'});
        ylabel('ADC');
        xlabel('Time (us)');
        grid on;

        linkaxes(ax_h, 'x');

        % Compress subplots upward to leave room for the global legend.
        for ch = 1:5
            pos    = get(ax_h(ch), 'Position');
            pos(2) = pos(2) * (1 - bot_margin) + bot_margin;
            pos(4) = pos(4) * (1 - bot_margin);
            set(ax_h(ch), 'Position', pos);
        end

        % Global horizontal legend in a hidden axes at the bottom.
        lgd_ax = axes('Parent', fig, ...
            'Position', [0.08, 0.005, 0.84, bot_margin - 0.01], ...
            'Visible', 'off');
        hold(lgd_ax, 'on');
        dummy_h    = gobjects(1, n_blocks);
        block_lbls = cell(1, n_blocks);
        for b = 1:n_blocks
            block_lbls{b} = sprintf('B%d', b - 1);
            dummy_h(b)    = plot(lgd_ax, NaN, NaN, '-', 'Color', cmap(b,:), 'LineWidth', 1.1);
        end
        legend(lgd_ax, dummy_h, block_lbls, ...
            'Orientation', 'horizontal', 'Location', 'north', 'Box', 'off');

        sgtitle(sprintf('Segment %d (%s)  rf->adc gap %.1fus, adc->adc gap %.1fus', ...
            s - 1, truth.base_name, seg.rf_adc_gap_us, seg.adc_adc_gap_us), 'Interpreter', 'none');
        figs(ii) = fig;
    end
end

function truth = coerce_truth(base_or_truth)
    if isstruct(base_or_truth)
        truth = base_or_truth;
    else
        truth = testutils.truth_parse_case(base_or_truth);
    end
end

    p = inputParser;
    addParameter(p, 'segment_idx', [], @(x) isempty(x) || isnumeric(x));
    parse(p, varargin{:});

    truth = coerce_truth(base_or_truth);

    GAMMA_HZ_PER_T = 42.577e6;           % proton gamma/(2pi), Hz/T
    GAMMA_GRAD     = GAMMA_HZ_PER_T * 1e-3;  % Hz/mT  -> mT/m = Hz/m  / GAMMA_GRAD
    GAMMA_RF       = GAMMA_HZ_PER_T * 1e-6;  % Hz/uT  ->  uT  = Hz    / GAMMA_RF

    if isempty(p.Results.segment_idx)
        seg_ids = 1:truth.segment_def.num_segments;
    else
        seg_ids = p.Results.segment_idx(:).';
    end

    figs = gobjects(1, numel(seg_ids));
    n_adc = double(truth.meta.adc_samples);

    for ii = 1:numel(seg_ids)
        s = seg_ids(ii);
        if s < 1 || s > truth.segment_def.num_segments
            error('truth:input', 'segment_idx %d out of range', s);
        end

        seg = truth.segment_def.segments(s);
        fig = figure('Name', sprintf('Segment %d: %s', s - 1, truth.base_name), 'Color', 'w');

        % Pre-compute shared block boundaries from max samples across all channels.
        block_width = zeros(1, seg.num_blocks);
        for b = 1:seg.num_blocks
            blk = seg.blocks(b);
            block_width(b) = max([double(blk.rf_n), double(blk.grad_n(:))']);
        end
        block_start = [0, cumsum(block_width + 1)];

        ax_h = gobjects(1, 5);
        grad_ylabels = {'Gx (mT/m)', 'Gy (mT/m)', 'Gz (mT/m)'};

        % --- Gradient channels ---
        for ch = 1:3
            ax_h(ch) = subplot(5, 1, ch);
            hold on;
            for b = 1:seg.num_blocks
                blk = seg.blocks(b);
                n = double(blk.grad_n(ch));
                if n > 0
                    y = double(blk.grad_amp(ch)) * double(blk.grad_wave{ch}(:)) / GAMMA_GRAD;
                    t = block_start(b) + double(blk.grad_delay(ch)) + (0:n-1).';
                    plot(t, y, 'LineWidth', 1.1, 'DisplayName', sprintf('B%d', b - 1));
                end
            end
            hold off;
            ylabel(grad_ylabels{ch});
            grid on;
            if ch == 1
                legend('show', 'Location', 'eastoutside');
            end
        end

        % --- RF |B1| ---
        ax_h(4) = subplot(5, 1, 4);
        hold on;
        for b = 1:seg.num_blocks
            blk = seg.blocks(b);
            n = double(blk.rf_n);
            if n > 0
                y = double(blk.rf_amp) * double(blk.rf_rho(:)) / GAMMA_RF;
                t = block_start(b) + double(blk.rf_delay) + (0:n-1).';
                plot(t, y, 'LineWidth', 1.1, 'DisplayName', sprintf('B%d', b - 1));
            end
        end
        hold off;
        ylabel('RF |B1| (uT)');
        grid on;

        % --- ADC mask ---
        ax_h(5) = subplot(5, 1, 5);
        hold on;
        for b = 1:seg.num_blocks
            blk = seg.blocks(b);
            if blk.has_adc && n_adc > 0
                t0 = block_start(b) + double(blk.adc_delay);
                t1 = t0 + n_adc - 1;
                fill([t0, t1, t1, t0], [0, 0, 1, 1], [0.2 0.6 0.9], ...
                    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end
        hold off;
        ylim([-0.1, 1.5]);
        yticks([0, 1]);
        yticklabels({'0', '1'});
        ylabel('ADC');
        xlabel('Pseudo-time (samples)');
        grid on;

        linkaxes(ax_h, 'x');

        sgtitle(sprintf('Segment %d (%s)  rf->adc gap %.3fus, adc->adc gap %.3fus', ...
            s - 1, truth.base_name, seg.rf_adc_gap_us, seg.adc_adc_gap_us), 'Interpreter', 'none');
        figs(ii) = fig;
    end
end

function truth = coerce_truth(base_or_truth)
    if isstruct(base_or_truth)
        truth = base_or_truth;
    else
        truth = testutils.truth_parse_case(base_or_truth);
    end
end
