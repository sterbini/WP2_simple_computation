# %%
import numpy as np

def beam_current_A (nb, Nb, frev):
    elementary_charge = 1.602176634e-19
    return nb*Nb*frev*elementary_charge

def sigma_x(normalized_emittance, beta, total_energy_GeV):
    proton_rest_mass_energy_GeV = 0.938272
    gamma_r = total_energy_GeV/proton_rest_mass_energy_GeV
    beta_r = np.sqrt(1 - 1/gamma_r**2)
    return np.sqrt(normalized_emittance*beta/gamma_r/beta_r)

def geometrical_emittance(normalized_emittance,  total_energy_GeV):
    proton_rest_mass_energy_GeV = 0.938272
    gamma_r = total_energy_GeV/proton_rest_mass_energy_GeV
    beta_r = np.sqrt(1 - 1/gamma_r**2)
    return normalized_emittance/gamma_r/beta_r

def geometrical_factor(sigma_x_m, full_crossing_angle_rad, sigma_z_m, full_crabbing_angle):
    theta_piwiski = (full_crossing_angle_rad+full_crabbing_angle)*sigma_z_m/2/sigma_x_m
    return 1/np.sqrt(1 + theta_piwiski**2)

def HO_tune_shift(N, normalized_emittance, geometrical_factor):
    # proton classical radius rp
    rp = 1.5347e-18
    # omitting the '-' sign
    return N*rp/normalized_emittance/4/np.pi*geometrical_factor

def beam_radius_IP(beta_star_m, geometrical_emittance):
    return np.sqrt(beta_star_m*geometrical_emittance)


def luminosity_loss_rate(bunch_intensity, number_of_bunches, luminosity, sigma_eff = 110e-31):
    return ((2 * luminosity**2 * sigma_eff) / (bunch_intensity*number_of_bunches))

def L_dLdt_hours(bunch_intensity, number_of_bunches, luminosity, sigma_eff = 110e-31):
    return luminosity/luminosity_loss_rate(bunch_intensity, number_of_bunches, luminosity, sigma_eff = sigma_eff)/3600.0

def print_it(my_string):
    print(6*'+', my_string, 6*'+')
    print(40*'*')
    print('Assuming the following parameters:')
    print(40*'*')

    print(f'Total Energy: {energy_total_Gev} GeV')
    print(f'Beta Star max: {beta_star_max*100} cm')
    print(f'Beta Star min: {beta_star_min*100} cm')
    print(f'Normalized Emittance: {emittance_normalized*1e6} um')
    print(f'Half Crossing Angle: {half_crossing_angle_rad*1e6} urad')
    print(f'Sigma Z: {sigma_z} m')
    print(f'Nb: {Nb/1e10}E10 ppb')
    print(f'nb: {nb} bunches')
    print(f'frequency of revolution: {frev} Hz')
    print(f'full crabbing angle: {full_crabbing_angle_rad*1e6} urad')

    print(40*'*')
    print('We get the following results: ')
    print(40*'*')

    print('Beam Current: ', beam_current_A(nb, Nb, frev)*1000, 'mA')
    print('Geometrical Emittance: ', geometrical_emittance(emittance_normalized, energy_total_Gev)*1e9, 'nm')
    print('Min Beam Radius at IP: ', beam_radius_IP(beta_star_min, geometrical_emittance(emittance_normalized, energy_total_Gev))*1e6, 'um')
    print('Geometrical Luminosity Factor: ', geometrical_factor_value)
    print(f'Head-On Tune Shift: {HO_tune_shift(Nb, emittance_normalized, geometrical_factor_value)*1e4}E-4')
    
    print(f'SOL L/dL/dt: {L_dldt_hours_SOL}')
    print(f'EOL L/dL/dt: {L_dldt_hours_EOL}')
    print(2*'\n')
# %% RunII 2015
energy_total_Gev = 6500
beta_star_max = 0.30
beta_star_min = 0.25
emittance_normalized = 2.1e-6
half_crossing_angle_rad = 160e-6
full_crossing_angle_rad = 2*half_crossing_angle_rad
sigma_z = 0.08
full_crabbing_angle_rad = 0
geometrical_factor_value = geometrical_factor(sigma_x(emittance_normalized, beta_star_max, energy_total_Gev), full_crossing_angle_rad, sigma_z, full_crabbing_angle_rad)
Nb = 1.1e11
nb = 2556
frev = 11245

nb_colliding_ip15 = 2544
bunch_intensity_SOL = 1.1e11
bunch_intensity_EOL = 0.9e11 # to be checked
luminosity = 2.1e38
L_dldt_hours_SOL = L_dLdt_hours(bunch_intensity_SOL, nb_colliding_ip15, luminosity, sigma_eff = 110e-31)
L_dldt_hours_EOL = L_dLdt_hours(bunch_intensity_EOL, nb_colliding_ip15, luminosity, sigma_eff = 110e-31)

print_it('Run2 2015 Parameters')

# %% RunIII 2022
energy_total_Gev = 6800
beta_star_max = 1.20
beta_star_min = 0.30
emittance_normalized = 2.1e-6
half_crossing_angle_rad = 160e-6
full_crossing_angle_rad = 2*half_crossing_angle_rad
sigma_z = 0.097
full_crabbing_angle_rad = 0
geometrical_factor_value = geometrical_factor(sigma_x(emittance_normalized, beta_star_max, energy_total_Gev), full_crossing_angle_rad, sigma_z, full_crabbing_angle_rad)
Nb = 1.6e11
nb = 2464
frev = 11245

nb_colliding_ip15 = 2452
bunch_intensity_SOL = 1.6e11
bunch_intensity_EOL = 1.28e11
luminosity = 2.0e38
L_dldt_hours_SOL = L_dLdt_hours(bunch_intensity_SOL, nb_colliding_ip15, luminosity, sigma_eff = 110e-31)
L_dldt_hours_EOL = L_dLdt_hours(bunch_intensity_EOL, nb_colliding_ip15, luminosity, sigma_eff = 110e-31)

print_it('Run3 2022 Parameters')

# %% HL-LHC
energy_total_Gev = 7000
beta_star_max = 1.00
beta_star_min = 0.15

emittance_normalized = 2.2e-6
half_crossing_angle_rad = 250e-6
full_crossing_angle_rad = 2*half_crossing_angle_rad
sigma_z = 0.076
full_crabbing_angle_rad = -190e-6*2
geometrical_factor_value = geometrical_factor(sigma_x(emittance_normalized, beta_star_max, energy_total_Gev),full_crossing_angle_rad, sigma_z, full_crabbing_angle_rad)
Nb = 2.2e11
nb = 2760
frev = 11245
crabbing_angle = -190e-6

nb_colliding_ip15 = 2748
bunch_intensity_SOL = 2.2e11
bunch_intensity_EOL = 1.22e11
luminosity = 5e38
L_dldt_hours_SOL = L_dLdt_hours(bunch_intensity_SOL, nb_colliding_ip15, luminosity, sigma_eff = 110e-31)
L_dldt_hours_EOL = L_dLdt_hours(bunch_intensity_EOL, nb_colliding_ip15, luminosity, sigma_eff = 110e-31)

print_it('HL-LHC Parameters')

# %%
