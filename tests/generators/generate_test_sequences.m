%% generate_test_sequences.m
%
% Iterative rebuild of segmentation test generation.
% Phase 1 scope:
%   - Build one basic GRE 2D case
%   - Export minimal segmentation-focused truth
%   - Keep placeholders for later TR/safety/freqmod truth

clear; clc;
import mr.*

% write_gre(true, 1, 1);
% write_gre(true, 3, 1);
% write_gre(true, 1, 3);
% write_gre(true, 3, 3);
% 
% write_mprage(true, 1, 1);
% write_mprage(true, 3, 1);
% write_mprage(true, 1, 3);
% write_mprage(true, 3, 3);

write_mprage_nav(true, 1, 1);
write_mprage_nav(true, 3, 1);
write_mprage_nav(true, 1, 3);
write_mprage_nav(true, 3, 3);

write_mprage_noncart(true, 1, 1, false);
write_mprage_noncart(true, 3, 1, false);
write_mprage_noncart(true, 1, 3, false);
write_mprage_noncart(true, 3, 3, false);


function seq = write_bssfp(write, num_slices, num_averages)
    base = sprintf('gre_2d_%dsl_%davg', num_slices, num_averages);
    fprintf('Generating sequence: %s\n', base);

    sys = make_system();

    % Basic GRE geometry intentionally small for fast iteration.
    fov = 0.22;
    Nx = 64;
    Ny = 8;
    slice_thickness = 5e-3;

    alpha = 10 * pi / 180;

    % --- create events ---
    [rf, gz, gz_reph] = mr.makeSincPulse(alpha * pi / 180, ...
        'Duration', rf_dur, 'SliceThickness', thick, ...
        'apodization', 0.5, 'timeBwProduct', 1.5, ...
        'system', sys, 'use', 'excitation');

    deltak = 1 / fov;
    gx = mr.makeTrapezoid('x', ...
        'FlatArea', Nx * deltak, ...
        'FlatTime', adc_dur, ...
        'system', sys);
    adc = mr.makeAdc(Nx, ...
        'Duration', gx.flatTime, ...
        'Delay', gx.riseTime, ...
        'system', sys);
    gx_pre = mr.makeTrapezoid('x', 'Area', -gx.area / 2, 'system', sys);
    phase_areas = ((0:Ny-1) - Ny/2) * deltak;

    % --- split & combine for bSSFP optimal timing ---
    gz_parts = mr.splitGradientAt(gz, mr.calcDuration(rf));
    gz_parts(1).delay = mr.calcDuration(gz_reph);
    gz_1 = mr.addGradients({gz_reph, gz_parts(1)}, 'system', sys);
    [rf, ~] = mr.align('right', rf, gz_1);
    gz_parts(2).delay = 0;
    gz_reph.delay = mr.calcDuration(gz_parts(2));
    gz_2 = mr.addGradients({gz_parts(2), gz_reph}, 'system', sys);

    gx_parts = mr.splitGradientAt(gx, ...
        ceil(mr.calcDuration(adc) / sys.gradRasterTime) * sys.gradRasterTime);
    gx_parts(1).delay = mr.calcDuration(gx_pre);
    gx_1 = mr.addGradients({gx_pre, gx_parts(1)}, 'system', sys);
    adc.delay = adc.delay + mr.calcDuration(gx_pre);
    gx_parts(2).delay = 0;
    gx_pre.delay = mr.calcDuration(gx_parts(2));
    gx_2 = mr.addGradients({gx_parts(2), gx_pre}, 'system', sys);

    pe_dur = mr.calcDuration(gx_2);

    gz_1.delay = max(mr.calcDuration(gx_2) - rf.delay + rf.ringdownTime, 0);
    rf.delay = rf.delay + gz_1.delay;

    TR = mr.calcDuration(gz_1) + mr.calcDuration(gx_1);

    % --- phase-encode template (max area, will be scaled) ---
    max_pe_area = max(abs(phase_areas));
    gyMax = mr.makeTrapezoid('y', 'Area', max_pe_area, 'Duration', pe_dur, 'system', sys);

    % --- pre-create labels ---
    lblOnce1 = mr.makeLabel('SET', 'ONCE', 1);
    lblOnce0 = mr.makeLabel('SET', 'ONCE', 0);
    lblOnce2 = mr.makeLabel('SET', 'ONCE', 2);

    % --- prep RF: half flip angle ---
    rf05 = rf;
    rf05.signal = 0.5 * rf.signal;

    % --- prep delay to center main acquisition around TR/2 ---
    prep_delay = mr.makeDelay( round((TR/2 - mr.calcDuration(gz_1)) / sys.gradRasterTime) * sys.gradRasterTime);
    gx_1_1 = mr.makeExtendedTrapezoidArea('x', 0, gx_2.first, -gx_2.area, sys);
    gy_pre_2 = mr.scaleGrad(gyMax, phase_areas(end) / max_pe_area);
    [prep_delay, gz_2, gy_pre_2, gx_1_1] = mr.align('left', prep_delay, gz_2, gy_pre_2, 'right', gx_1_1);

    for z = 1:num_slices
        rf05.freqOffset = gz.amplitude * thick * (z - 1 - (num_slices-1)/2);
        rf.freqOffset = gz.amplitude * thick * (z - 1 - (num_slices-1)/2);
        
        % --- alpha/2 prep (ONCE=1) ---
        seq.addBlock(rf05, gz_1, lblOnce1);
        seq.addBlock(prep_delay, gz_2, gy_pre_2, gx_1_1);

        % --- main loop ---
        for i = 1:Ny
            rf.phaseOffset = pi * mod(i, 2);
            adc.phaseOffset = pi * mod(i, 2);

            gy_pre_1 = mr.scaleGrad(gy_pre_2, -1);             % undo previous PE
            gy_pre_2 = mr.scaleGrad(gyMax, phase_areas(i) / max_pe_area);  % new PE

            if i == 1
                seq.addBlock(rf, gz_1, gy_pre_1, gx_2, lblOnce0); % clear ONCE flag -> first main block
            else
                seq.addBlock(rf, gz_1, gy_pre_1, gx_2);
            end
            seq.addBlock(gx_1, gy_pre_2, gz_2, adc);
        end

        % --- exit block (ONCE=2) ---
        seq.addBlock(gx_2, lblOnce2);
    end


    seq.setDefinition('FOV', [fov fov slice_thickness * num_slices]);
    seq.setDefinition('NumSlices', num_slices);

    [ok, err] = seq.checkTiming;
    if ~ok
        error('Timing check failed:\n%s', strjoin(err, '\n'));
    end

    if write == false
        return;
    end

    out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');

    tb = TruthBuilder(seq, sys);
    tb.setBlocksPerTR(3 + 2*Ny);
    tb.setSegments([3 + 2*Ny]);
    tb.setSegmentOrder([1]);
    tb.setNumAverages(num_averages);
    tb.export(out_dir, base);
end


function seq = write_gre(write, num_slices, num_averages)
    base = sprintf('gre_2d_%dsl_%davg', num_slices, num_averages);
    fprintf('Generating sequence: %s\n', base);

    sys = make_system();

    % Basic GRE geometry intentionally small for fast iteration.
    fov = 0.22;
    Nx = 64;
    Ny = 8;
    slice_thickness = 5e-3;
    ndummy = 5;

    alpha = 10 * pi / 180;
    rf_spoil_inc = 84.0; % degrees

    % RF and slice-select
    [rf, gz] = mr.makeSincPulse(alpha, ...
        'Duration', 2.0e-3, ...
        'SliceThickness', slice_thickness, ...
        'timeBwProduct', 4, ...
        'apodization', 0.5, ...
        'use', 'excitation', ...
        'system', sys);
    gz_reph = mr.makeTrapezoid('z', 'Area', -gz.area/2, 'Duration', 1.0e-3, 'system', sys);
    gz_spoil = mr.makeTrapezoid('z', 'Area', 4 / slice_thickness, 'Duration', 1.0e-3, 'system', sys);

    % Readout and ADC
    readout_time = 2.56e-3;
    gx_full = mr.makeTrapezoid('x', 'FlatArea', Nx/fov, 'FlatTime', readout_time, 'system', sys);
    gx_parts = mr.splitGradientAt(gx_full, gx_full.riseTime + gx_full.flatTime);
    gx = gx_parts(1); % truncate at end of flat
    adc = mr.makeAdc(Nx, 'Duration', gx_full.flatTime, 'Delay', gx_full.riseTime, 'system', sys);
    dummy_adc = mr.makeDelay(mr.calcDuration(adc));

    % Pre/rewinder templates
    gx_pre = mr.makeTrapezoid('x', 'Area', -gx_full.area/2, 'Duration', 1.0e-3, 'system', sys);
    gx_spoil_area = 4 / slice_thickness;

    % Bridged spoiler that starts at gx flat amplitude for continuity.
    gx_spoil = mr.makeExtendedTrapezoidArea('x', gx_full.amplitude, 0, gx_spoil_area, sys);

    pe_areas = ((0:Ny-1) - floor(Ny/2)) / fov;
    max_pe_area = max(abs(pe_areas));
    gy_phase = mr.makeTrapezoid('y', 'Area', max_pe_area, 'Duration', 1.0e-3, 'system', sys);

    % ONCE labels for prep/main semantics.
    lbl_once1 = mr.makeLabel('SET', 'ONCE', 1);
    lbl_once0 = mr.makeLabel('SET', 'ONCE', 0);

    seq = mr.Sequence(sys);

    rf_center = mr.calcRfCenter(rf);
    rf_phase = 0.0;
    rf_inc = 0.0;

    % Dummy TRs (once region): same 4-block TR shape with PE=0.
    for d = 1:ndummy
        rf_inc = mod(rf_inc + rf_spoil_inc, 360.0);
        rf_phase = mod(rf_phase + rf_inc, 360.0);

        rf_curr = rf;
        rf_curr.phaseOffset = rf_phase / 180 * pi;

        gy_pre = mr.scaleGrad(gy_phase, 0.0);
        gy_rew = mr.scaleGrad(gy_phase, 0.0);

        if d == 1
            seq.addBlock(rf_curr, gz, lbl_once1);
        else
            seq.addBlock(rf_curr, gz);
        end
        seq.addBlock(gx_pre, gy_pre, gz_reph);
        seq.addBlock(gx, dummy_adc);
        seq.addBlock(gx_spoil, gy_rew, gz_spoil);
    end

    % Main imaging loop: PE -> slices (slice is inner loop).
    first_main = true;
    for pe = 1:Ny
        if max_pe_area > 0
            yscale = pe_areas(pe) / max_pe_area;
        else
            yscale = 0;
        end
        for sl = 1:num_slices
            rf_inc = mod(rf_inc + rf_spoil_inc, 360.0);
            rf_phase = mod(rf_phase + rf_inc, 360.0);

            rf_curr = rf;
            slc_shift = (sl - 1 - (num_slices - 1) / 2);

            rf_curr.freqOffset = gz.amplitude * slice_thickness * slc_shift;
            rf_curr.phaseOffset = rf_phase / 180 * pi - 2 * pi * rf_curr.freqOffset * rf_center;

            adc_curr = adc;
            adc_curr.freqOffset = rf_curr.freqOffset;
            adc_curr.phaseOffset = rf_phase / 180 * pi;

            gy_pre = mr.scaleGrad(gy_phase, yscale);
            gy_rew = mr.scaleGrad(gy_phase, -yscale);

            if first_main
                seq.addBlock(rf_curr, gz, lbl_once0);
                first_main = false;
            else
                seq.addBlock(rf_curr, gz);
            end
            seq.addBlock(gx_pre, gy_pre, gz_reph);
            seq.addBlock(gx, adc_curr);
            seq.addBlock(gx_spoil, gy_rew, gz_spoil);
        end
    end

    seq.setDefinition('FOV', [fov fov slice_thickness * num_slices]);
    seq.setDefinition('NumSlices', num_slices);

    [ok, err] = seq.checkTiming;
    if ~ok
        error('Timing check failed:\n%s', strjoin(err, '\n'));
    end

    if write == false
        return;
    end

    out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');

    tb = TruthBuilder(seq, sys);
    tb.setBlocksPerTR(4);
    tb.setSegments([4]);
    tb.setSegmentOrder([1]);
    tb.setNumAverages(num_averages);
    tb.export(out_dir, base);
end


function seq = write_mprage(write, num_slices, num_averages)
    base = sprintf('mprage_2d_%dsl_%davg', num_slices, num_averages);
    fprintf('Generating sequence: %s\n', base);

    sys = make_system();

    % Basic GRE geometry intentionally small for fast iteration.
    fov = 0.22;
    Nx = 64;
    Ny = 8;
    slice_thickness = 5e-3;

    alpha = 10 * pi / 180;
    rf_spoil_inc = 84.0; % degrees
    
    % Inversion
    % RF and slice-select
    rf180 = mr.makeBlockPulse(pi, ...
        'Duration', 20.0e-3, ...
        'use', 'excitation', ... % should be 'inversion', but this way we get timing
        'system', sys);
    delayTI = mr.makeDelay(0.1e-3);

    % RF and slice-select
    [rf, gz] = mr.makeSincPulse(alpha, ...
        'Duration', 2.0e-3, ...
        'SliceThickness', slice_thickness, ...
        'timeBwProduct', 4, ...
        'apodization', 0.5, ...
        'use', 'excitation', ...
        'system', sys);
    gz_reph = mr.makeTrapezoid('z', 'Area', -gz.area/2, 'Duration', 1.0e-3, 'system', sys);
    gz_spoil = mr.makeTrapezoid('z', 'Area', 4 / slice_thickness, 'Duration', 1.0e-3, 'system', sys);

    % Readout and ADC
    readout_time = 2.56e-3;
    gx_full = mr.makeTrapezoid('x', 'FlatArea', Nx/fov, 'FlatTime', readout_time, 'system', sys);
    gx_parts = mr.splitGradientAt(gx_full, gx_full.riseTime + gx_full.flatTime);
    gx = gx_parts(1); % truncate at end of flat
    adc = mr.makeAdc(Nx, 'Duration', gx_full.flatTime, 'Delay', gx_full.riseTime, 'system', sys);

    % Pre/rewinder templates
    gx_pre = mr.makeTrapezoid('x', 'Area', -gx_full.area/2, 'Duration', 1.0e-3, 'system', sys);
    gx_spoil_area = 4 / slice_thickness;

    % Bridged spoiler that starts at gx flat amplitude for continuity.
    gx_spoil = mr.makeExtendedTrapezoidArea('x', gx_full.amplitude, 0, gx_spoil_area, sys);
    delayTR = mr.makeDelay(0.5e-3);

    pe_areas = ((0:Ny-1) - floor(Ny/2)) / fov;
    max_pe_area = max(abs(pe_areas));
    gy_phase = mr.makeTrapezoid('y', 'Area', max_pe_area, 'Duration', 1.0e-3, 'system', sys);

    seq = mr.Sequence(sys);

    rf_center = mr.calcRfCenter(rf);
    rf_phase = 0.0;
    rf_inc = 0.0;

    % Main imaging loop: PE -> slices (slice is inner loop).
    for sl = 1:num_slices
        slc_shift = (sl - 1 - (num_slices - 1) / 2);
        
        seq.addBlock(rf180);
        seq.addBlock(gz_spoil);
        seq.addBlock(delayTI);
        
        for pe = 1:Ny
            if max_pe_area > 0
                yscale = pe_areas(pe) / max_pe_area;
            else
                yscale = 0;
            end
        
            rf_inc = mod(rf_inc + rf_spoil_inc, 360.0);
            rf_phase = mod(rf_phase + rf_inc, 360.0);

            rf_curr = rf;
            rf_curr.freqOffset = gz.amplitude * slice_thickness * slc_shift;
            rf_curr.phaseOffset = rf_phase / 180 * pi - 2 * pi * rf_curr.freqOffset * rf_center;

            adc_curr = adc;
            adc_curr.freqOffset = rf_curr.freqOffset;
            adc_curr.phaseOffset = rf_phase / 180 * pi;

            gy_pre = mr.scaleGrad(gy_phase, yscale);
            gy_rew = mr.scaleGrad(gy_phase, -yscale);

            seq.addBlock(rf_curr, gz);
            seq.addBlock(gx_pre, gy_pre, gz_reph);
            seq.addBlock(gx, adc_curr);
            seq.addBlock(gx_spoil, gy_rew, gz_spoil);
        end
        seq.addBlock(delayTR);
    end

    seq.setDefinition('FOV', [fov fov slice_thickness * num_slices]);
    seq.setDefinition('NumSlices', num_slices);

    [ok, err] = seq.checkTiming;
    if ~ok
        error('Timing check failed:\n%s', strjoin(err, '\n'));
    end

    if write == false
        return;
    end

    out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');

    tb = TruthBuilder(seq, sys);
    tb.setBlocksPerTR(2 + 1 + 4 * Ny + 1);
    tb.setSegments([2, 1, 4]);
    tb.setSegmentOrder([1, 2, 3 * ones(1, Ny), 2]);
    tb.setNumAverages(num_averages);
    tb.export(out_dir, base);
end


function seq = write_mprage_nav(write, num_slices, num_averages)
    base = sprintf('mprage_nav_2d_%dsl_%davg', num_slices, num_averages);
    fprintf('Generating sequence: %s\n', base);

    sys = make_system();

    % Basic GRE geometry intentionally small for fast iteration.
    fov = 0.22;
    Nx = 64;
    Ny = 8;
    slice_thickness = 5e-3;

    alpha = 10 * pi / 180;
    rf_spoil_inc = 84.0; % degrees
    
    % Inversion
    % RF and slice-select
    rf180 = mr.makeBlockPulse(pi, ...
        'Duration', 20.0e-3, ...
        'use', 'excitation', ... % should be 'inversion', but this way we get timing
        'system', sys);
    delayTI = mr.makeDelay(0.1e-3);

    % RF and slice-select
    [rf, gz] = mr.makeSincPulse(alpha, ...
        'Duration', 2.0e-3, ...
        'SliceThickness', slice_thickness, ...
        'timeBwProduct', 4, ...
        'apodization', 0.5, ...
        'use', 'excitation', ...
        'system', sys);
    gz_reph = mr.makeTrapezoid('z', 'Area', -gz.area/2, 'Duration', 1.0e-3, 'system', sys);
    gz_spoil = mr.makeTrapezoid('z', 'Area', 4 / slice_thickness, 'Duration', 1.0e-3, 'system', sys);

    % Readout and ADC
    readout_time = 2.56e-3;
    gx_full = mr.makeTrapezoid('x', 'FlatArea', Nx/fov, 'FlatTime', readout_time, 'system', sys);
    gx_parts = mr.splitGradientAt(gx_full, gx_full.riseTime + gx_full.flatTime);
    gx = gx_parts(1); % truncate at end of flat
    adc = mr.makeAdc(Nx, 'Duration', gx_full.flatTime, 'Delay', gx_full.riseTime, 'system', sys);

    % Navigator: 3 orthogonal single-shot spirals (16 samples)
    % Each nav has a slice-selective sinc pulse selecting the respective plane.
    Nx_nav = 16;
    nav_slice_thickness = 10e-3;
    [gx_nav_cells, gy_nav_cells, adc_nav] = testutils.makeTestSpiral(sys, 1, Nx_nav, 10.0*fov);
    
    gx_nav_ax = gx_nav_cells{1}; gx_nav_rew_ax = mr.makeExtendedTrapezoidArea('x', gx_nav_ax.last, 0.0, -gx_nav_ax.area, sys);
    gy_nav_ax = gy_nav_cells{1}; gy_nav_rew_ax = mr.makeExtendedTrapezoidArea('y', gy_nav_ax.last, 0.0, -gy_nav_ax.area, sys);
    gz_nav_spoil_ax = mr.makeTrapezoid('z', 'Area', 4 / nav_slice_thickness, 'Duration', 1.0e-3, 'system', sys);

    gx_nav_cor = gx_nav_ax; gx_nav_cor.channel = 'x'; gx_nav_rew_cor = mr.makeExtendedTrapezoidArea('x', gx_nav_cor.last, 0.0, -gx_nav_cor.area, sys);
    gz_nav_cor = gy_nav_ax; gz_nav_cor.channel = 'z'; gz_nav_rew_cor = mr.makeExtendedTrapezoidArea('z', gz_nav_cor.last, 0.0, -gz_nav_cor.area, sys);
    gy_nav_spoil_cor = mr.makeTrapezoid('y', 'Area', 4 / nav_slice_thickness, 'Duration', 1.0e-3, 'system', sys);
    
    gy_nav_sag = gx_nav_ax; gy_nav_sag.channel = 'y'; gy_nav_rew_sag = mr.makeExtendedTrapezoidArea('y', gy_nav_sag.last, 0.0, -gy_nav_sag.area, sys);
    gz_nav_sag = gy_nav_ax; gz_nav_sag.channel = 'z'; gz_nav_rew_sag = mr.makeExtendedTrapezoidArea('z', gz_nav_sag.last, 0.0, -gz_nav_sag.area, sys);
    gx_nav_spoil_sag = mr.makeTrapezoid('x', 'Area', 4 / nav_slice_thickness, 'Duration', 1.0e-3, 'system', sys);

    % Axial nav: slice-select on z
    [rf_nav_ax, gz_nav_ss] = mr.makeSincPulse(alpha, ...
        'Duration', 0.5e-3, ...
        'SliceThickness', nav_slice_thickness, ...
        'timeBwProduct', 4, ...
        'apodization', 0.5, ...
        'use', 'excitation', ...
        'system', sys);
    
    % Coronal nav: slice-select on y
    gy_nav_ss = gz_nav_ss; gy_nav_ss.channel = 'y';
    rf_nav_cor = rf_nav_ax;
    
    % Sagittal nav: slice-select on x
    gx_nav_ss = gz_nav_ss; gx_nav_ss.channel = 'x';
    rf_nav_sag = rf_nav_ax;

    lblNav = mr.makeLabel('SET', 'NAV', 1);

    % Pre/rewinder templates
    gx_pre = mr.makeTrapezoid('x', 'Area', -gx_full.area/2, 'Duration', 1.0e-3, 'system', sys);
    gx_spoil_area = 4 / slice_thickness;

    % Bridged spoiler that starts at gx flat amplitude for continuity.
    gx_spoil = mr.makeExtendedTrapezoidArea('x', gx_full.amplitude, 0, gx_spoil_area, sys);
    delayTR = mr.makeDelay(0.5e-3);

    pe_areas = ((0:Ny-1) - floor(Ny/2)) / fov;
    max_pe_area = max(abs(pe_areas));
    gy_phase = mr.makeTrapezoid('y', 'Area', max_pe_area, 'Duration', 1.0e-3, 'system', sys);

    seq = mr.Sequence(sys);

    rf_center = mr.calcRfCenter(rf);
    rf_phase = 0.0;
    rf_inc = 0.0;

    % Main imaging loop: PE -> slices (slice is inner loop).
    for sl = 1:num_slices
        slc_shift = (sl - 1 - (num_slices - 1) / 2);
        
        seq.addBlock(rf180);
        seq.addBlock(gz_spoil);
        seq.addBlock(delayTI);
        
        for pe = 1:Ny
            if max_pe_area > 0
                yscale = pe_areas(pe) / max_pe_area;
            else
                yscale = 0;
            end
        
            rf_inc = mod(rf_inc + rf_spoil_inc, 360.0);
            rf_phase = mod(rf_phase + rf_inc, 360.0);

            rf_curr = rf;
            rf_curr.freqOffset = gz.amplitude * slice_thickness * slc_shift;
            rf_curr.phaseOffset = rf_phase / 180 * pi - 2 * pi * rf_curr.freqOffset * rf_center;

            adc_curr = adc;
            adc_curr.freqOffset = rf_curr.freqOffset;
            adc_curr.phaseOffset = rf_phase / 180 * pi;

            gy_pre = mr.scaleGrad(gy_phase, yscale);
            gy_rew = mr.scaleGrad(gy_phase, -yscale);

            seq.addBlock(rf_curr, gz);
            seq.addBlock(gx_pre, gy_pre, gz_reph);
            seq.addBlock(gx, adc_curr);
            seq.addBlock(gx_spoil, gy_rew, gz_spoil);
        end
        % Navigator: 3 orthogonal single-shot spirals
        % Axial (xy plane): slice-select on z
        seq.addBlock(rf_nav_ax, gz_nav_ss, lblNav);
        seq.addBlock(gx_nav_ax, gy_nav_ax, adc_nav);
        seq.addBlock(gx_nav_rew_ax, gy_nav_rew_ax, gz_nav_spoil_ax);
        
        % Coronal (xz plane): slice-select on y
        seq.addBlock(rf_nav_cor, gy_nav_ss, lblNav);
        seq.addBlock(gx_nav_cor, gz_nav_cor, adc_nav);
        seq.addBlock(gx_nav_rew_cor, gy_nav_spoil_cor, gz_nav_rew_cor);
        
        % Sagittal (yz plane): slice-select on x
        seq.addBlock(rf_nav_sag, gx_nav_ss, lblNav);
        seq.addBlock(gy_nav_sag, gz_nav_sag, adc_nav);
        seq.addBlock(gx_nav_rew_ax, gy_nav_rew_ax);
        seq.addBlock(gx_nav_spoil_sag, gy_nav_rew_sag, gz_nav_rew_sag);
 
        seq.addBlock(delayTR);
    end

    seq.setDefinition('FOV', [fov fov slice_thickness * num_slices]);
    seq.setDefinition('NumSlices', num_slices);

    [ok, err] = seq.checkTiming;
    if ~ok
        error('Timing check failed:\n%s', strjoin(err, '\n'));
    end

    if write == false
        return;
    end

    out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');

    tb = TruthBuilder(seq, sys);
    tb.setBlocksPerTR(2 + 1 + 4 * Ny + 3 * 2 + 1);
    tb.setSegments([2, 1, 4, 2, 2, 2]);
    tb.setSegmentOrder([1, 2, 3 * ones(1, Ny), 4, 5, 6, 2]);
    tb.setNumAverages(num_averages);
    tb.export(out_dir, base);
end


function seq = write_mprage_noncart(write, Nz, num_averages, use_rotext)
    base = sprintf('mprage_noncart_3d_%dsl_%davg_userotext%d', Nz, num_averages, use_rotext);
    fprintf('Generating sequence: %s\n', base);

    sys = make_system();

    % Basic geometry intentionally small for fast iteration.
    fov = 0.22;
    Nx = 64;
    num_shots = 8;
    slab_thickness = 15e-3;

    alpha = 10 * pi / 180;
    rf_spoil_inc = 84.0; % degrees

    % Inversion
    rf180 = mr.makeBlockPulse(pi, ...
        'Duration', 20.0e-3, ...
        'use', 'excitation', ...
        'system', sys);
    delayTI = mr.makeDelay(0.1e-3);

    % Excitation
    rf = mr.makeBlockPulse(alpha, ...
        'Duration', 2.0e-3, ...
        'use', 'excitation', ...
        'system', sys);
    gz_spoil = mr.makeTrapezoid('z', 'Area', 4 / 5.0e-3, 'Duration', 1.0e-3, 'system', sys);

    % --- Spiral readout (multi-interleave) ---
    if use_rotext
        % Rotext path: generate only first interleave; rotation events
        % will rotate it to the other interleave angles.
        [gx_cells, gy_cells, adc] = testutils.makeTestSpiral(sys, 1, Nx, fov);
        gx_base = gx_cells{1};
        gy_base = gy_cells{1};
        dphi = 2 * pi / num_shots;  % interleave angular step
    else
        % Explicit path: pre-compute all interleaves across all partitions.
        [gx_shots, gy_shots, adc] = testutils.makeTestSpiral(sys, num_shots * Nz, Nx, fov);
        spoke_idx = 0;
    end

    % Partition encoding (along z)
    gz_phase = mr.makeTrapezoid('z', 'Area', -Nz / fov / 2, 'system', sys);

    delayTR = mr.makeDelay(0.5e-3);

    if Nz > 1
        z_areas = ((0:Nz-1) - floor(Nz/2)) / slab_thickness;
        max_z_area = max(abs(z_areas));
    else
        z_areas = 0;
        max_z_area = 0;
    end

    seq = mr.Sequence(sys);

    rf_center = mr.calcRfCenter(rf);
    rf_phase = 0.0;
    rf_inc = 0.0;

    % Main imaging loop: partitions -> shots
    for z = 1:Nz
        if max_z_area > 0
            zscale = z_areas(z) / max_z_area;
        else
            zscale = 0;
        end

        seq.addBlock(rf180);
        seq.addBlock(gz_spoil);
        seq.addBlock(delayTI);

        for i = 1:num_shots
            rf_inc = mod(rf_inc + rf_spoil_inc, 360.0);
            rf_phase = mod(rf_phase + rf_inc, 360.0);

            rf_curr = rf;
            rf_curr.phaseOffset = rf_phase / 180 * pi - 2 * pi * rf_curr.freqOffset * rf_center;

            adc_curr = adc;
            adc_curr.phaseOffset = rf_phase / 180 * pi;

            seq.addBlock(rf_curr);
            if use_rotext
                rot_angle = (i - 1) * dphi;
                seq.addBlock(adc_curr, gx_base, gy_base, ...
                    mr.scaleGrad(gz_phase, zscale), ...
                    mr.makeRotation('axis', 'z', 'angle', rot_angle));
            else
                spoke_idx = spoke_idx + 1;
                seq.addBlock(adc_curr, gx_shots{spoke_idx}, gy_shots{spoke_idx}, ...
                    mr.scaleGrad(gz_phase, zscale));
            end
            seq.addBlock(gz_spoil);
        end

        seq.addBlock(delayTR);
    end

    seq.setDefinition('FOV', [fov fov slab_thickness]);
    seq.setDefinition('NumPartitions', Nz);

    [ok, err] = seq.checkTiming;
    if ~ok
        error('Timing check failed:\n%s', strjoin(err, '\n'));
    end

    if write == false
        return;
    end

    out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');

    tb = TruthBuilder(seq, sys);
    tb.setBlocksPerTR(2 + 1 + 3 * num_shots + 1);
    tb.setSegments([2, 1, 3]);
    tb.setSegmentOrder([1, 2, 3 * ones(1, num_shots), 2]);
    tb.setNumAverages(num_averages);
    tb.export(out_dir, base);
end


function sys = make_system()
    sys = mr.opts( ...
        'MaxGrad',   28,   'GradUnit', 'mT/m', ...
        'MaxSlew',   150,  'SlewUnit', 'T/m/s', ...
        'rfRasterTime',         2e-6, ...
        'gradRasterTime',      20e-6, ...
        'adcRasterTime',        2e-6, ...
        'blockDurationRaster', 20e-6);
end
