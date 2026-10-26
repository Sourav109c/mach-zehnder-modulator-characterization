"""Synthetic MZM characterization: 8-panel dashboard (not experimental data)."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import lfilter

OUT = Path(__file__).resolve().parent / 'results'
OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(42)
WAVELENGTH_NM = 1550.0
PIN_MW = 1.0
TRUE_VPI = 4.0
TRUE_IL_DB = 3.0
TRUE_ER_DB = 25.0
BIAS_V = 0.35
BANDWIDTH_GHZ = 3.0  # assumed small-signal electrical-to-optical 3-dB bandwidth
BIT_RATE_GBPS = 1.0
SAMPLES_PER_BIT = 200

def mzm_power(voltage, vpi, il_db, er_db, bias_v):
    """Balanced intensity MZM, single-arm voltage convention."""
    minimum_fraction = 10.0 ** (-er_db / 10.0)
    phase = np.pi * (np.asarray(voltage) - bias_v) / vpi
    return PIN_MW * 10.0 ** (-il_db / 10.0) * (
        minimum_fraction + (1 - minimum_fraction) * np.cos(phase / 2) ** 2
    )

voltage = np.linspace(0, 10, 501)
ideal = PIN_MW * np.cos(np.pi * voltage / (2 * TRUE_VPI)) ** 2
actual = mzm_power(voltage, TRUE_VPI, TRUE_IL_DB, TRUE_ER_DB, BIAS_V)
v_meas = np.linspace(0, 10, 161)
p_meas = np.clip(mzm_power(v_meas, TRUE_VPI, TRUE_IL_DB, TRUE_ER_DB, BIAS_V)
                 + RNG.normal(0, 0.0012, len(v_meas)), 1e-9, None)
fit, _ = curve_fit(mzm_power, v_meas, p_meas,
                   p0=[3.8, 3.2, 20, 0.2],
                   bounds=([2, 0, 10, -1], [6, 10, 40, 1]), maxfev=30000)
fit_vpi, fit_il, fit_er, fit_bias = fit
fit_curve = mzm_power(voltage, *fit)
residuals = p_meas - mzm_power(v_meas, *fit)
rmse = np.sqrt(np.mean(residuals ** 2))

# OOK: drive between constructive and destructive interference voltages.
bit_pattern = np.array([1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0])
bit_duration_ns = 1 / BIT_RATE_GBPS
time_ns = np.arange(len(bit_pattern) * SAMPLES_PER_BIT) * bit_duration_ns / SAMPLES_PER_BIT
bit_index = np.repeat(np.arange(len(bit_pattern)), SAMPLES_PER_BIT)
drive_v = fit_bias + fit_vpi * bit_pattern[bit_index]
ideal_ook = mzm_power(drive_v, *fit)

# First-order low-pass response H(f) = 1/(1+j*f/f3dB).
# Filtering the voltage before the nonlinear transfer function approximates
# the small-signal EO response; it is NOT a full traveling-wave RF model.
dt_ns = bit_duration_ns / SAMPLES_PER_BIT
tau_ns = 1 / (2 * np.pi * BANDWIDTH_GHZ)
a = np.exp(-dt_ns / tau_ns)
# Initialize at the first drive voltage to avoid artificial startup transient.
filtered_drive = lfilter([1 - a], [1, -a], drive_v - drive_v[0]) + drive_v[0]
limited_ook = mzm_power(filtered_drive, *fit)
noisy_ook = np.clip(limited_ook + RNG.normal(0, 0.004, len(limited_ook)), 0, None)

freq_ghz = np.linspace(0, 10, 500)
s21_db = -10 * np.log10(1 + (freq_ghz / BANDWIDTH_GHZ) ** 2)

fig, axs = plt.subplots(2, 4, figsize=(21, 10), constrained_layout=True)
fig.suptitle('Mach–Zehnder Modulator | Synthetic Characterization Dashboard',
             fontsize=18, fontweight='bold')

ax = axs[0, 0]
ax.plot(voltage, ideal, label='Ideal (lossless)')
ax.plot(voltage, actual, label='Finite IL and ER')
ax.set(title='1. Optical transfer function', xlabel='Voltage (V)', ylabel='Power (mW)')
ax.legend(fontsize=8)

ax = axs[0, 1]
ax.scatter(v_meas, p_meas, s=10, alpha=0.6, label='Synthetic measurements')
ax.plot(voltage, fit_curve, label='Fitted model')
ax.set(title='2. Model fitting', xlabel='Voltage (V)', ylabel='Power (mW)')
ax.legend(fontsize=8)

ax = axs[0, 2]
ax.scatter(v_meas, residuals * 1000, s=10)
ax.axhline(0, color='gray', linestyle='--')
ax.set(title=f'3. Fit residuals | RMSE = {rmse*1000:.2f} µW',
       xlabel='Voltage (V)', ylabel='Residual (µW)')

ax = axs[0, 3]
ax.plot(voltage, np.pi * (voltage - fit_bias) / fit_vpi)
ax.axhline(np.pi, linestyle='--', color='gray', label='π rad')
ax.set(title='4. Differential phase shift', xlabel='Voltage (V)', ylabel='Phase (rad)')
ax.legend(fontsize=8)

ax = axs[1, 0]
ax.step(time_ns, drive_v, where='post')
ax.set(title='5. Electrical OOK drive (1 Gb/s)', xlabel='Time (ns)',
       ylabel='Voltage (V)', xlim=(0, 16))

ax = axs[1, 1]
ax.plot(time_ns, ideal_ook, label='Instantaneous')
ax.plot(time_ns, limited_ook, '--', linewidth=1.5, label='3 GHz low-pass')
ax.set(title='6. Optical OOK output', xlabel='Time (ns)',
       ylabel='Power (mW)', xlim=(0, 16))
ax.legend(fontsize=8)

ax = axs[1, 2]
ax.plot(freq_ghz, s21_db)
ax.axhline(-3.0103, color='gray', linestyle='--', label='−3.01 dB')
ax.axvline(BANDWIDTH_GHZ, color='tab:red', linestyle=':', label=f'{BANDWIDTH_GHZ:g} GHz')
ax.set(title='7. EO S21 magnitude (model)', xlabel='Frequency (GHz)',
       ylabel='Normalized |S21| (dB)', xlim=(0, 10))
ax.legend(fontsize=8)

ax = axs[1, 3]
# Eye: fold waveform into 2-bit windows; discard the first two bits to
# avoid initialization artifacts. Start windows on bit boundaries.
for start_bit in range(2, len(bit_pattern) - 1):
    start = start_bit * SAMPLES_PER_BIT
    segment = noisy_ook[start:start + 2 * SAMPLES_PER_BIT]
    if len(segment) == 2 * SAMPLES_PER_BIT:
        ax.plot(np.arange(len(segment)) / SAMPLES_PER_BIT, segment,
                color='tab:blue', alpha=0.27, linewidth=0.8)
ax.set(title='8. Optical eye diagram (3 GHz + noise)',
       xlabel='Time (UI)', ylabel='Power (mW)', xlim=(0, 2))

for ax in axs.flat:
    ax.grid(True, alpha=0.25)

summary = (
    'MZM CHARACTERIZATION — SYNTHETIC DATA ONLY\n'
    f'Wavelength: {WAVELENGTH_NM:.0f} nm\n'
    f'Input optical power: {PIN_MW:.3f} mW\n'
    f'Fitted Vpi: {fit_vpi:.3f} V (true {TRUE_VPI:.3f} V)\n'
    f'Fitted insertion loss: {fit_il:.3f} dB (true {TRUE_IL_DB:.3f} dB)\n'
    f'Fitted extinction ratio: {fit_er:.3f} dB (true {TRUE_ER_DB:.3f} dB)\n'
    f'Fitted bias offset: {fit_bias:.3f} V (true {BIAS_V:.3f} V)\n'
    f'Fit RMSE: {rmse * 1000:.3f} µW\n'
    f'Assumed small-signal EO 3-dB bandwidth: {BANDWIDTH_GHZ:.2f} GHz\n'
    f'OOK bit rate: {BIT_RATE_GBPS:.2f} Gb/s\n'
    'Notes: EO S21 is an assumed first-order low-pass model, not a measured VNA trace.\n'
    'Eye diagram uses a short illustrative bit sequence with simulated noise.\n'
    'No experimental hardware data or measured bandwidth are used.\n'
)
print('\n' + summary)
fig.savefig(OUT / 'mzm_dashboard_8graphs.png', dpi=170, bbox_inches='tight')
np.savetxt(OUT / 'synthetic_voltage_sweep.csv', np.column_stack((v_meas, p_meas)),
           delimiter=',', header='voltage_V,power_mW', comments='')
np.savetxt(OUT / 'synthetic_s21.csv', np.column_stack((freq_ghz, s21_db)),
           delimiter=',', header='frequency_GHz,S21_dB', comments='')
(OUT / 'characterization_summary.txt').write_text(summary, encoding='utf-8')
print(f'Results saved to: {OUT}')
plt.show()
