function figs = truth_plot_segments(base_or_truth, varargin)
%TRUTH_PLOT_SEGMENTS Plot segment block waveforms on pseudo-time.
%
% Segment defs do not store explicit sample times for every waveform. This
% plot uses pseudo-time (sample index plus delay offset) to inspect shape,
% continuity, and per-block amplitude scaling.

    p = inputParser;
    addParameter(p, 'segment_idx', [], @(x) isempty(x) || isnumeric(x));
    parse(p, varargin{:});

    truth = coerce_truth(base_or_truth);

    if isempty(p.Results.segment_idx)
        seg_ids = 1:truth.segment_def.num_segments;
    else
        seg_ids = p.Results.segment_idx(:).';
    end

    figs = gobjects(1, numel(seg_ids));

    for ii = 1:numel(seg_ids)
        s = seg_ids(ii);
        if s < 1 || s > truth.segment_def.num_segments
            error('truth:input', 'segment_idx %d out of range', s);
        end

        seg = truth.segment_def.segments(s);
        fig = figure('Name', sprintf('Segment %d: %s', s - 1, truth.base_name), 'Color', 'w');

        chan_names = {'Gx', 'Gy', 'Gz', 'RF'};
        chan_idx = [1, 2, 3];

        for ax = 1:3
            subplot(4, 1, ax);
            hold on;
            cursor = 0;
            for b = 1:seg.num_blocks
                blk = seg.blocks(b);
                n = double(blk.grad_n(chan_idx(ax)));
                if n > 0
                    y = double(blk.grad_amp(chan_idx(ax))) * double(blk.grad_wave{chan_idx(ax)}(:));
                    t = cursor + double(blk.grad_delay(chan_idx(ax))) + (0:n-1).';
                    plot(t, y, 'LineWidth', 1.1, 'DisplayName', sprintf('B%d', b - 1));
                end
                cursor = cursor + max(1, n + 1);
            end
            hold off;
            ylabel(chan_names{ax});
            grid on;
            if ax == 1
                legend('show', 'Location', 'eastoutside');
            end
        end

        subplot(4, 1, 4);
        hold on;
        cursor = 0;
        for b = 1:seg.num_blocks
            blk = seg.blocks(b);
            n = double(blk.rf_n);
            if n > 0
                y = double(blk.rf_amp) * double(blk.rf_rho(:));
                t = cursor + double(blk.rf_delay) + (0:n-1).';
                plot(t, y, 'LineWidth', 1.1, 'DisplayName', sprintf('B%d', b - 1));
            end
            cursor = cursor + max(1, n + 1);
        end
        hold off;
        ylabel('RF |B1|');
        xlabel('Pseudo-time (samples)');
        grid on;

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
