from species import SpeciesInit
from species.data.database import Database
from species.fit.dretrieval import AtmosphericRetrieval
from species.plot.plot_mcmc import plot_posterior
from species.plot.plot_retrieval import plot_clouds, plot_opacities, plot_pt_profile
from species.plot.plot_spectrum import plot_spectrum
from species.util.fit_util import get_residuals

SpeciesInit()
database = Database()

# import urllib.request
# urllib.request.urlretrieve('http://irtfweb.ifa.hawaii.edu/~spex/IRTF_Spectral_Library/Data/L3_2MASSJ1506+1321.fits',
#                            'L3_2MASSJ1506+1321.fits')

# database.add_object('2MASS J15065441+1321060',
#                     parallax=(85.4250, 0.1902),
#                     app_mag=None,
#                     flux_density=None,
#                     spectrum={'IRTF': ('L3_2MASSJ1506+1321.fits', None, 2000.)},
#                     deredden=None)

database.add_object('AF Lep b',
                    parallax=(37.2539,0.0195),
                    app_mag=None,
                    flux_density=None,
                    spectrum={'GRAVITY': ('balmer2024b_aflepb_gravity_kband_contrast-absflux-andcovars_final.fits', 'balmer2024b_aflepb_gravity_kband_contrast-absflux-andcovars_final.fits', 500.)},
                    deredden=None)


retrieve = AtmosphericRetrieval(object_name='AF Lep b',
                                line_species=['CO_all_iso_HITEMP', 'H2O_HITEMP', 'CH4', 'NH3', 'CO2', 'Na_allard', 'K_allard', 'TiO_all_Exomol', 'VO_Plez', 'FeH', 'H2S'],
                                cloud_species=['MgSiO3(c)_cd'],
                                output_folder='multinest',
                                wavel_range=(1.9, 2.5),
                                scattering=True,
                                inc_spec=True,
                                inc_phot=False,
                                pressure_grid='smaller',
                                weights=None,
                                ccf_species=None,
                                max_pressure=1e3,
                                lbl_opacity_sampling=None)

retrieve.run_dynesty(bounds={'logg': (2.5, 6.0),
                               'c_o_ratio': (0.1, 1.5),
                               'metallicity': (-0.3, 3.),
                               'radius': (1.1, 2.),
                               'fsed': (1., 8.),
                               'log_kzz': (8., 12.),
                               'sigma_lnorm': (1.2, 5.)},
                       chemistry='equilibrium',
                       quenching='pressure',
                       pt_profile='eddington',
                       fit_corr=None,
                       n_live_points=100,
                       npool=16,
                       dynamic=False,
                       mpi_pool=False,
                       dlogz=0.5,
                       resume=True,
                       plotting=True,
                       pt_smooth=0.3,
                       temp_nodes=None,
                       abund_nodes=None,
                       prior={'c_o_ratio': (0.6, 0.05),
                              'metallicity': (1., 0.01),
                              'radius': (1.3, 0.05),
                              'mass': (3.5, 0.5),
                              'log_kzz': (10., 1.),
                              'log_p_quench': (2.3, 1.),
                              },
                       check_phot_press=None)

database.add_retrieval(tag='aflep_gravonly',
                       output_folder='multinest',
                       inc_teff=False)

_, _ = database.get_retrieval_teff(tag='aflep_gravonly',
                                   random=30)

fig = plot_posterior(tag='aflep_gravonly',
                             offset=(-0.3, -0.35),
                             vmr=False,
                             inc_luminosity=False,
                             inc_mass=False,
                             inc_pt_param=False,
                             inc_loglike=False,
                             output='aflep_gravonly_post.pdf')

samples, radtrans = database.get_retrieval_spectra(tag='aflep_gravonly',
                                                   random=30,
                                                   wavel_range=(0.5, 5.),
                                                   spec_res=500.)

fig = plot_pt_profile(tag='aflep_gravonly',
                              random=100,
                              xlim=(500., 6000.),
                              offset=(-0.07, -0.14),
                              output='aflep_gravonly_ptprofile.pdf',
                              radtrans=radtrans,
                              extra_axis='grains',
                              rad_conv_bound=False)

fig = plot_opacities(tag='aflep_gravonly',
                    offset=(-0.1, -0.14),
                    output='aflep_gravonly_opacities.pdf',
                    radtrans=radtrans)

object_box = database.get_object('AF Lep b')

best = database.get_median_sample(tag='aflep_gravonly')

model_box = radtrans.get_model(model_param=best,
                               spec_res=500.,
                               wavel_resample=None,
                               plot_contribution=False)

res_box = get_residuals(datatype='model',
                        spectrum='petitradtrans',
                        parameters=best,
                        objectbox=object_box,
                        inc_phot=True,
                        inc_spec=True,
                        radtrans=radtrans)

fig = plot_spectrum(boxes=[samples, model_box, object_box],
                            filters=None,
                            plot_kwargs=[{'ls': '-', 'lw': 0.1, 'color': 'gray'},
                                         {'ls': '-', 'lw': 0.5, 'color': 'black'},
                                         {'GRAVITY': {'ls': '-', 'lw': 0.3, 'color': 'dodgerblue', 'label': 'GRAVITY'}}],
                            residuals=res_box,
                            # xlim=(0.8, 4.2),
                            # ylim=(0., 2.7e-14),
                            # ylim_res=(-22., 22.),
                            scale=('linear', 'linear'),
                            offset=(-0.6, -0.05),
                            figsize=(10, 6),
                            legend=[{'loc': 'upper right', 'fontsize': 10.}, None],
                            output='aflep_gravonly_spectrum.pdf')
