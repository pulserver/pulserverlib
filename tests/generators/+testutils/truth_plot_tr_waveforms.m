function fig = truth_plot_tr_waveforms(base_or_truth, varargin)
%TRUTH_PLOT_TR_WAVEFORMS Plot canonical TR waveforms.
%
%   truth_plot_tr_waveforms('gre_2d_1sl_1avg')
%   truth_plot_tr_waveforms(truth_struct, 'overlay', true)

    p = inputParser;
    addParameter(p, 'overlay', true, @islogical);
    addParameter(p, 'show_slew', false, @islogical);
    parse(p, varargin{:});

    truth = coerce_truth(base_or_truth);

    fig = figure('Name', sprintf('Canonical TR waveforms: %s', truth.base_name), 'Color', 'w');
    tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    labels = {'Gx', 'Gy', 'Gz'};
    fields = {'gx', 'gy', 'gz'};

    for ax = 1:3
        nexttile;
        hold on;
        for i = 1:truth.tr_waveforms.num_trs
            w = truth.tr_waveforms.waveforms(i);
            t_ms = double(w.time_us) / 1e3;
            y = double(w.(fields{ax}));
            if p.Results.show_slew
                dt = diff(t_ms) * 1e-3;
                dy = diff(y);
                y = [0; dy ./ max(dt, eps)];
            end
            if p.Results.overlay
                plot(t_ms, y, 'LineWidth', 1.2, 'DisplayName', sprintf('TR %d', i - 1));
            else
                plot(t_ms, y, 'LineWidth', 1.2);
            end
        end
        hold off;
        ylabel(labels{ax});
        grid on;
        if ax == 1 && p.Results.overlay
            legend('show', 'Location', 'best');
        end
        if ax == 3
            xlabel('Time (ms)');
        end
    end

    title(tl, sprintf('Canonical TRs (%s)', truth.base_name), 'Interpreter', 'none');
end

function truth = coerce_truth(base_or_truth)
    if isstruct(base_or_truth)
        truth = base_or_truth;
    else
        truth = truth_parse_case(base_or_truth);
    end
end
