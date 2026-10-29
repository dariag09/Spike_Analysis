function Spike_analysis_help
% SPIKE_ANALYSIS_HELP  User guide for Spike_analysis.m (patch-clamp AP analysis).
%
% PURPOSE
%   Spike_analysis.m is a MATLAB script for detecting and characterizing
%   action potentials (APs) in voltage recordings stored in Axon Binary
%   Files (.abf). It analyzes one selected channel across all sweeps,
%   estimates AP waveform features, and plots per-spike measurements.
%   Original script attribution: Daria Garkun.
%
% QUICK START
%   1. Add Spike_analysis.m and the external abfload function to the
%      MATLAB path. Ensure the ABF file contains voltage recordings in mV.
%   2. Open Spike_analysis.m and adjust the constants near the top.
%   3. Run Spike_analysis in MATLAB; select an .abf file in the dialog.
%   4. Review the Command Window sweep counts and generated figures.
%   5. Inspect workspace variables for the final sweep; the script does
%      not currently write a spreadsheet or save results automatically.
%
% INPUTS AND DEPENDENCIES
%   Input: Axon .abf file selected interactively with uigetfile.
%   Loader: [data,si,h] = abfload(dataFileName) (external dependency).
%   Expected data dimensions: samples x channels x sweeps.
%   si / h.si: sampling interval in microseconds.
%   h.dataPtsPerChan: samples per channel in each sweep.
%   Voltage values are assumed to be in mV; sample intervals in us.
%   MATLAB version must support scripts with local functions and tiledlayout.
%   There is no cancellation handling in the file-selection dialog.
%
% USER-EDITABLE PARAMETERS (values in supplied script)
%   AP_Criteria          = 0     mV     Minimum accepted AP peak voltage.
%                                       Also used to locate positive-voltage
%                                       episodes for the preliminary ISI pass.
%   AP_Threshold         = 20    mV/ms  Minimum positive dV/dt for AP onset.
%                                       Numerically equal to V/s.
%   DataChannel          = 3     index  Channel to analyze (1-based).
%   AP_TimeToPeack       = 1.5   ms     Search interval after detected onset
%                                       for AP maximum (spelling as in code).
%   window_slope_search  = 20    ms     Upper bound for post-onset/post-peak
%                                       search windows, shortened according
%                                       to minimum preliminary interspike
%                                       interval across the file.
%   window_before_sweep  = 4     ms     Converted to sample indices but
%                                       otherwise unused in current script.
%   window_after_sweep   = 4     ms     Converted to sample indices but
%                                       otherwise unused in current script.
%
% PROCESSING PIPELINE
%   A. Load the ABF file and construct a time vector in seconds.
%   B. For each sweep, find distinct runs with voltage > AP_Criteria;
%      findspike returns the maximum-voltage sample of each run.
%      Compute the shortest difference between successive detected peak
%      sample indices (preliminary interspike interval, in samples).
%   C. Set idx_slope_search to the smaller of the configured slope
%      search duration and the shortest preliminary interspike interval.
%   D. For each sweep, compute the first-difference voltage derivative:
%          derivative_Sweep = [diff(Sweep); 0] * (1000 / h.si)
%      with Sweep in mV and h.si in microseconds. This yields mV/ms.
%   E. findonset identifies the first sample in each contiguous run
%      where derivative_Sweep > AP_Threshold. These are candidate onsets.
%   F. Search the next AP_TimeToPeack milliseconds for the voltage peak;
%      accept a candidate only when its peak exceeds AP_Criteria.
%   G. Compute the waveform metrics below for each accepted event.
%   H. Create a figure for sweeps with more than two derivative-detected
%      candidate events (the figure condition does not count only accepted
%      events).
%
% MEASUREMENTS (for each accepted AP)
%   AP_Threshold_Idx      Sample index of derivative-based onset.
%   AP_Threshold_Time     Onset time (s).
%   AP_Threshold_Amp      Voltage at onset (mV); despite 'Amp' name,
%                         this is an absolute membrane potential.
%   AP_Peak_Idx           Index of maximum voltage in peak-search window.
%   AP_Peak_Time          Time of that maximum (s).
%   AP_Peak_Voltage       Voltage at the peak (mV); stored in a matrix
%                         indexed by spike and sweep.
%   AP_Peak_Amp           Intended AP height = peak voltage - onset
%                         voltage (mV). See implementation warning below.
%   AP_rise_slope         Largest dV/dt from onset through slope window
%                         (mV/ms, numerically V/s).
%   AP_decay_slope        Most negative dV/dt from peak through slope
%                         window (mV/ms, numerically V/s).
%   AP_Ahp_Voltage        Minimum voltage after peak within slope window
%                         (mV).
%   AP_Ahp_Amp            AP_Ahp_Voltage - AP_Threshold_Amp (mV);
%                         typically negative if AHP is below onset.
%   AP_Halfamp_Duration   Time between first and last samples above
%                         half-height in the post-onset search interval
%                         (ms). This is sample-based, not interpolated.
%   interspikeInterval    Minimum preliminary peak-to-peak separation
%                         per sweep (samples), not converted to seconds.
%   interspikeInterval_Min  Minimum valid interval across sweeps (samples).
%
% FORMULAS
%   AP height (intended) = Vpeak - Vonset.
%   Half-height voltage = Vonset + (Vpeak - Vonset)/2.
%   Half-width (ms) = (last_above_half - first_above_half) * si/1000.
%   AHP relative to onset = Vminimum_postpeak - Vonset.
%   Maximum rise/decay slopes are extrema of the sampled derivative,
%   not fitted slopes.
%
% FIGURE CONTENT
%   Each qualifying sweep generates a 3-by-6 tiled figure with:
%     - Membrane voltage versus time, green onset and red peak markers.
%     - Voltage-derivative phase-plane plot (dV/dt versus voltage),
%       labeled 'Polar plot' in the original code.
%     - Per-AP traces of onset voltage, AP height, maximum rise slope,
%       maximum decay slope, AHP relative to onset, and half-width.
%   The voltage trace is displayed with fixed x limits [0.4, 2.5] s
%   and y limits [-100, 50] mV; half-width panel limits [0, 0.4] ms.
%
% OUTPUT AND STORAGE
%   The script has no formal output arguments and does not save figures
%   or tabular results. Many AP measurement arrays are reset on each
%   sweep, so the final workspace does not retain complete per-sweep
%   results. AP_Peak_Voltage is indexed by sweep but is not cleared in
%   the per-sweep reset. The script calls clear all and close all at start.
%
% IMPLEMENTATION WARNINGS / REVIEW BEFORE QUANTITATIVE USE
%   1. AP_Peak_Amp uses AP_Peak_Voltage(end,1), rather than the current
%      spike/sweep's peak. For multiple spikes or sweeps this can give
%      incorrect AP heights and, consequently, half-height/half-width.
%   2. findspike and findonset use sparse-index bookkeeping that should
%      be tested on synthetic traces; multi-event boundaries may be
%      incorrectly detected. findspike returns NaN when no event exists.
%   3. The preliminary interval computation uses size(event_Idx,2),
%      assumes a row-vector return, and may mishandle missing events.
%      If no sweep contains two spikes, the global search interval can
%      become empty, making later indexing fail.
%   4. Candidate peak/slope windows are not clipped to sweep length;
%      events near the end of a sweep may trigger indexing errors.
%   5. Millisecond-to-index conversions are not explicitly rounded;
%      fractional sample counts can cause indexing errors.
%   6. Half-width takes the first and last above-half-height samples
%      in a search window. If none exist, indexing fails; multiple
%      crossings or overlapping spikes can distort the measurement.
%   7. Search windows are limited by the shortest preliminary ISI over
%      all sweeps, potentially shortening AHP and decay measurements.
%   8. A plot is requested after >2 candidate onsets even if fewer
%      than three pass the peak criterion. Plotting may then be empty.
%   9. AHP is measured relative to onset, not to resting membrane
%      potential; its sign differs from some positive-magnitude usages.
%  10. This code has been documented from static inspection; the
%      measurement calculations and plotting have not been validated
%      against example ABF recordings.
%
% LOCAL FUNCTIONS
%   findspike(Wave, EventThreshold)
%       Finds above-voltage-threshold runs and returns their peak indices.
%       Used for preliminary spike separation, not the final AP onset.
%   findonset(Wave, EventThreshold)
%       Finds onset indices of above-threshold derivative runs.
%       Used to propose AP onset locations.
%
% EXAMPLE
%   >> help Spike_analysis_help
%   >> Spike_analysis
%
% SEE ALSO
%   Spike_analysis, abfload, diff, max, min, tiledlayout.
%
% NOTE
%   This companion file provides help text only; it does not execute
%   or modify Spike_analysis.m. Use 'help Spike_analysis_help' to read it.
end
