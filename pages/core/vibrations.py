'''
                    Waveplot: An online wavefunction viewer
                    Copyright (C) 2023  Jon G. C. Kragskow

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
import numpy as np
from numpy.typing import NDArray, ArrayLike
from functools import lru_cache
from dash import dcc, html, Input, Output, State, callback, no_update, Patch
import dash_bootstrap_components as dbc
from . import common as com
from scipy.special import factorial
import scipy.constants as con
import io
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import copy


# c in cm s-1
LIGHT = con.speed_of_light * 100
HBAR = con.Planck

VIB_LAYOUT = copy.copy(com.BASIC_LAYOUT)
VIB_LAYOUT.xaxis.title = {
    'text': 'x (Å)',
    'font': {
        'family': 'Arial',
        'size': 14,
        'color': 'black'
    }
}
VIB_LAYOUT.yaxis.title = {
    'text': 'Energy (cm<sup>-1</sup>)',
    'font': {
        'family': 'Arial',
        'size': 14,
        'color': 'black'
    }
}

VIB_LAYOUT.showlegend = False
VIB_LAYOUT.yaxis.tickformat = '.2f'


def hermite(n: int, x: ArrayLike) -> NDArray:
    '''
    Computes nth hermite polynomial values for a given set of x values

    Parameters
    ----------
    n: int
        Order of polynomial
    x: array_like
        Values of x used in H^n (x)

    Returns
    -------
    ndarray of floats
        Hermite polynomial of nth order evaluated at x
    '''

    if n == 0:
        return x * 0 + 1
    elif n == 1:
        return 2 * x
    else:
        return 2 * x * hermite(n - 1, x) - 2 * (n-1) * hermite(n - 2, x)


def calc_harmonic_energies(k: float, m: float, max_n: int) -> tuple[
                               NDArray, NDArray, NDArray, float]:
    '''
    Calculate energy of harmonic oscillator as both classical and quantum
    entity

    Parameters
    ----------
    k : float
        Force constant (N/m)
    m : float
        Reduced mass (kg)
    max_n : int
        maximum value of n used for Harmonic states
    Returns
    -------
    ndarray of floats
        Harmonic state energies for quantum oscillator in Joules
    ndarray of floats
        Harmonic energies for classical oscillator in Joules
    ndarray of floats
        Displacements used for classical oscillator in metres
    float
        Zero point displacement in metres
    '''

    # Angular frequency
    omega = np.sqrt(k / m) / (2*np.pi)  # rad s^-1

    # Harmonic state energies
    state_E = np.array([HBAR * omega * (n + 0.5) for n in range(0, max_n + 1)])

    # Harmonic potential
    # E = 1/2 kx^2
    max_x = np.sqrt((max_n + 0.5) * 2 * HBAR * omega / k)

    displacement = np.linspace(-max_x*1.2, max_x*1.2, 100)  # m
    harmonic_E = 0.5 * k * displacement**2  # J

    # Find zero point displacement
    zpd = np.sqrt(HBAR * omega / k)  # m

    return state_E, harmonic_E, displacement, zpd


def calculate_mu(w, k):
    '''
    Calculates 'reduced mass' from angular frequency and force constant

    Parameters
    ----------
    w : float
        Angular frequency Omega (s^-1)
    k : float
        Force constant k (N m^-1)

    Returns
    -------
    float
        Reduced mass mu (g mol^-1)
    '''

    mu = k / w**2
    # Convert mass to kg
    mu /= 1.6605E-27

    return mu


def calculate_k(w, mu):
    '''
    Calculates force constant from angular frequency and reduced mass

    Parameters
    ----------
    w : float
        Angular frequency Omega (s^-1)
    mu : float
        Reduced mass (g mol^-1)

    Returns
    -------
    float
        Force constant k (N m^-1)
    '''

    # Convert mass to kg
    mu *= 1.6605E-27
    k = mu * w**2

    return k


def harmonic_wf(n: int, x: ArrayLike, m: float, omega: float):
    '''
    Calculates normalised harmonic wavefunction for nth state over x

    Parameters
    ----------
    n : int
        Harmonic oscillator quantum number
    x : array_like
        Displacement (m)
    m : float
        Reduced mass (kg)
    omega : float
        Angular frequency (float)
    Returns
    -------
    ndarray of floats
        Harmonic wavefunction of nth levels evaluated at each point in x
    '''

    x = np.asarray(x)

    beta = np.sqrt(m*omega / HBAR)

    # Compute hermite polynomial
    h = hermite(n, beta * x)

    # Normalisation factor
    N = 1. / np.sqrt(2**n * factorial(n)) * ((m * omega)/(np.pi * HBAR))**0.25

    # Wavefunction
    wf = h * N * np.exp(-beta**2 * x**2 * 0.5)
    print(n, h)

    return wf


class OptionsDiv(com.Div):
    def __init__(self, prefix, **kwargs):
        # Initialise base class attributes
        super().__init__(prefix=prefix, **kwargs)

        self.lin_wn_input = dbc.Input(
            id=self.prefix('lin_wn'),
            placeholder=2888,
            type='number',
            min=0.0001,
            value=2888,
            style={
                'textAlign': 'center'
            }
        )
        self.lin_wn_check = dbc.Checkbox(
            value=False,
            id=self.prefix('lin_wn_fix')
        )
        self.lin_wn_ig = self.make_input_group(
            [
                dbc.InputGroupText(u'v\u0334'),
                self.lin_wn_input,
                dbc.InputGroupText(r'cm⁻¹'),
                dbc.InputGroupText(
                    self.lin_wn_check
                )
            ]
        )

        self.ang_wn_input = dbc.Input(
            id=self.prefix('ang_wn'),
            placeholder=18145.84,
            type='number',
            min=0.0001,
            value=18145.84,
            style={
                'textAlign': 'center'
            }
        )
        self.ang_wn_check = dbc.Checkbox(
            value=False,
            id=self.prefix('ang_wn_fix')
        )
        self.ang_wn_ig = self.make_input_group(
            [
                dbc.InputGroupText(u'ῶ'),
                self.ang_wn_input,
                dbc.InputGroupText(r'cm⁻¹'),
                dbc.InputGroupText(
                    self.ang_wn_check
                )
            ]
        )

        self.fc_input = dbc.Input(
            id=self.prefix('fc'),
            placeholder=480,
            type='number',
            min=0.,
            value=480,
            style={
                'textAlign': 'center'
            }
        )
        self.fc_check = dbc.Checkbox(
            value=True,
            id=self.prefix('fc_fix')
        )
        self.fc_ig = self.make_input_group(
            [
                dbc.InputGroupText(u'k'),
                self.fc_input,
                dbc.InputGroupText(r'N m⁻¹'),
                dbc.InputGroupText(
                    self.fc_check
                )
            ]
        )

        self.rm_input = dbc.Input(
            id=self.prefix('rm'),
            placeholder=0.9768,
            type='number',
            min=0.,
            value=0.9768,
            style={
                'textAlign': 'center'
            }
        )
        self.rm_check = dbc.Checkbox(
            value=True,
            id=self.prefix('rm_fix')
        )
        self.rm_ig = self.make_input_group(
            [
                dbc.InputGroupText(u'μ'),
                self.rm_input,
                dbc.InputGroupText(r'g mol⁻¹'),
                dbc.InputGroupText(
                    self.rm_check
                )
            ]
        )

        self.var_store = dcc.Store(
            id=self.prefix('var_store'),
            data={}
        )

        self.max_n_input = dbc.Input(
            id=self.prefix('max_n'),
            placeholder=5,
            value=5,
            type='number',
            min=0,
            max=25,
            style={
                'textAlign': 'center'
            }
        )

        self.max_n_ig = self.make_input_group(
            [
                dbc.InputGroupText('Max. n'),
                self.max_n_input
            ]
        )

        self.wf_scale_input = dbc.Input(
            id=self.prefix('wf_scale'),
            placeholder=5,
            value=5,
            type='number',
            min=0,
            max=25,
            style={
                'textAlign': 'center'
            }
        )

        self.wf_scale_ig = self.make_input_group(
            [
                dbc.InputGroupText('WF Scale'),
                self.wf_scale_input
            ]
        )

        self.text_size_input = dbc.Input(
            id=self.prefix('text_size'),
            placeholder=18,
            type='number',
            min=10,
            max=25,
            value=15,
            style={
                'textAlign': 'center',
                'verticalAlign': 'middle',
                'horizontalAlign': 'middle'
            }
        )

        self.wf_colour_select = dbc.Select(
            id=self.prefix('wf_colour'),
            options=[
                {
                    'label': 'Default',
                    'value': 'normal'
                },
                {
                    'label': 'Tol',
                    'value': 'tol'
                },
                {
                    'label': 'Wong',
                    'value': 'wong'
                }
            ],
            value='normal',
            style={
                'textAlign': 'center',
                'verticalAlign': 'middle',
                'horizontalAlign': 'middle',
                'alignItems': 'auto',
                'display': 'inline'
            }
        )

        self.wf_colour_ig = self.make_input_group(
            [
                dbc.InputGroupText('WF colour'),
                self.wf_colour_select
            ]
        )

        self.wf_linewidth_input = dbc.Input(
            id=self.prefix('wf_linewidth'),
            placeholder=3,
            type='number',
            min=1,
            max=10,
            value=3,
            style={
                'textAlign': 'center',
                'verticalAlign': 'middle',
                'horizontalAlign': 'middle'
            }
        )

        self.wf_linewidth_ig = self.make_input_group(
            [
                dbc.InputGroupText('Linewidth'),
                self.wf_linewidth_input
            ]
        )

        self.tog_wf_check = dbc.Checkbox(
            value=True,
            id=self.prefix('toggle_wf')
        )
        self.wf_toggle_ig = self.make_input_group(
            [
                dbc.InputGroupText('Toggle'),
                dbc.InputGroupText(
                    self.tog_wf_check
                )
            ]
        )

        self.wf_modal = self.make_modal(
            [
                self.wf_colour_ig,
                self.wf_linewidth_ig,
                self.wf_toggle_ig
            ]
        )

        self.wf_modal_btn = dbc.Button(
            'WF options',
            id=self.prefix('open_wf_modal'),
            color='primary',
            className='me-1',
            n_clicks=0,
            style={
                'textAlign': 'center',
                'width': '100%'
            }
        )

        self.pe_colour_select = dbc.Input(
            id=self.prefix('PE_colour'),
            type='color',
            value='#ffffff',
            style={
                'height': '40px'
            }
        )

        self.pe_colour_ig = self.make_input_group(
            [
                dbc.InputGroupText('pe colour'),
                self.pe_colour_select
            ]
        )

        self.pe_linewidth_input = dbc.Input(
            id=self.prefix('PE_linewidth'),
            placeholder=3,
            type='number',
            min=1,
            max=10,
            value=3,
            style={
                'textAlign': 'center',
                'verticalAlign': 'middle',
                'horizontalAlign': 'middle'
            }
        )

        self.pe_linewidth_ig = self.make_input_group(
            [
                dbc.InputGroupText('Linewidth'),
                self.pe_linewidth_input
            ]
        )

        self.tog_pe_check = dbc.Checkbox(
            value=True,
            id=self.prefix('toggle_pe')
        )
        self.pe_toggle_ig = self.make_input_group(
            [
                dbc.InputGroupText('Toggle'),
                dbc.InputGroupText(
                    self.tog_pe_check
                )
            ]
        )

        self.pe_modal = self.make_modal(
            [
                self.pe_colour_ig,
                self.pe_linewidth_ig,
                self.pe_toggle_ig
            ]
        )

        self.pe_modal_btn = dbc.Button(
            'PE options',
            id=self.prefix('open_pe_modal'),
            color='primary',
            className='me-1',
            n_clicks=0,
            style={
                'textAlign': 'center',
                'width': '100%'
            }
        )

        self.state_colour_select = dbc.Input(
            id=self.prefix('state_colour'),
            type='color',
            value='#ffffff',
            style={
                'height': '40px'
            }
        )

        self.state_colour_ig = self.make_input_group(
            [
                dbc.InputGroupText('state colour'),
                self.state_colour_select
            ]
        )

        self.state_linewidth_input = dbc.Input(
            id=self.prefix('state_linewidth'),
            placeholder=3,
            type='number',
            min=1,
            max=10,
            value=3,
            style={
                'textAlign': 'center',
                'verticalAlign': 'middle',
                'horizontalAlign': 'middle'
            }
        )

        self.state_linewidth_ig = self.make_input_group(
            [
                dbc.InputGroupText('Linewidth'),
                self.state_linewidth_input
            ]
        )

        self.tog_state_check = dbc.Checkbox(
            value=True,
            id=self.prefix('toggle_state')
        )
        self.state_toggle_ig = self.make_input_group(
            [
                dbc.InputGroupText('Toggle'),
                dbc.InputGroupText(
                    self.tog_state_check
                )
            ]
        )

        self.state_modal = self.make_modal(
            [
                self.state_colour_ig,
                self.state_linewidth_ig,
                self.state_toggle_ig
            ]
        )

        self.state_modal_btn = dbc.Button(
            'State options',
            id=self.prefix('open_state_modal'),
            color='primary',
            className='me-1',
            n_clicks=0,
            style={
                'textAlign': 'center',
                'width': '100%'
            }
        )

        self.text_size_ig = self.make_input_group(
            [
                dbc.InputGroupText('Text size'),
                self.text_size_input
            ]
        )

        self.download_button = dbc.Button(
            'Download Data',
            id=self.prefix('download_data'),
            style={
                'boxShadow': 'none',
                'textalign': 'top'
            }
        )
        self.download_trigger = dcc.Download(
            id=self.prefix('download_data_trigger')
        )

        self.image_format_select = dbc.Select(
            id=self.prefix('save_format'),
            style={
                'textAlign': 'center',
                'horizontalAlign': 'center',
                'display': 'inline'
            },
            options=[
                {
                    'label': 'svg',
                    'value': 'svg',
                },
                {
                    'label': 'png',
                    'value': 'png',
                },
                {
                    'label': 'jpeg',
                    'value': 'jpeg',
                }
            ],
            value='svg'
        )

        self.image_format_ig = self.make_input_group(
            [
                dbc.InputGroupText('Image format'),
                self.image_format_select
            ]
        )

        self.make_div_contents()
        return

    def make_input_group(self, elements):

        group = dbc.InputGroup(
            elements,
            class_name='mb-3',
        )
        return group

    def make_modal(self, elements):

        modal = dbc.Modal(
            elements,
            class_name='mb-3',
        )
        return modal

    def make_div_contents(self):
        '''
        Assembles div children in rows and columns
        '''

        contents = [
            dbc.Row([
                dbc.Col(
                    html.H4(
                        id=self.prefix('parameters_header'),
                        style={
                            'textAlign': 'center',
                        },
                        children='Parameters'
                    )
                ),
                dbc.Tooltip(
                    children='Toggle buttons specify the two fixed parameters',
                    target=self.prefix('parameters_header'),
                    style={
                        'textAlign': 'center',
                    }
                )
            ]),
            dbc.Row(
                [
                    dbc.Col(
                        self.lin_wn_ig,
                        class_name='mb-3'
                    ),
                    dbc.Col(
                        self.ang_wn_ig,
                        class_name='mb-3'
                    )
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        self.fc_ig,
                        class_name='mb-3'
                    ),
                    dbc.Col(
                        self.rm_ig,
                        class_name='mb-3'
                    ), self.var_store
                ]
            ),
            dbc.Row([
                dbc.Col(
                    html.H4(
                        style={
                            'textAlign': 'center',
                        },
                        children='Plot options'
                    )
                )
            ]),
            dbc.Row([
                dbc.Col(
                    self.max_n_ig
                )
            ]),
            dbc.Row([
                dbc.Col(
                    [
                        self.wf_modal_btn,
                        self.wf_modal
                    ],
                    class_name='mb-3'
                ),
                dbc.Col(
                    [
                        self.pe_modal_btn,
                        self.pe_modal
                    ],
                    class_name='mb-3'
                ),
                dbc.Col(
                    [
                        self.state_modal_btn,
                        self.state_modal
                    ],
                    class_name='mb-3'
                ),
            ]),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            self.download_button,
                            self.download_trigger
                        ],
                        class_name='mb-3'
                    ),
                    dbc.Col(
                        self.image_format_ig,
                        class_name='mb-3'
                    )
                ],
                className='align-items-center'
            )
        ]

        self.div.children = contents
        return


def assemble_callbacks(plot_div: com.PlotDiv, options: OptionsDiv):

    # Toggle editable inputs
    outputs = [
        Output(options.var_store, 'data'),
        Output(options.lin_wn_input, 'value'),
        Output(options.ang_wn_input, 'value'),
        Output(options.fc_input, 'value'),
        Output(options.rm_input, 'value'),
        Output(options.lin_wn_input, 'disabled'),
        Output(options.ang_wn_input, 'disabled'),
        Output(options.fc_input, 'disabled'),
        Output(options.rm_input, 'disabled')
    ]
    inputs = [
        Input(options.lin_wn_input, 'value'),
        Input(options.ang_wn_input, 'value'),
        Input(options.fc_input, 'value'),
        Input(options.rm_input, 'value'),
        Input(options.lin_wn_check, 'value'),
        Input(options.ang_wn_check, 'value'),
        Input(options.fc_check, 'value'),
        Input(options.rm_check, 'value')
    ]

    callback(outputs, inputs)(update_inputs)

    # Calculate data and store
    inputs = [
        Input(options.var_store, 'data'),
        Input(options.max_n_input, 'value')
    ]
    callback(
        Output(plot_div.store, 'data'),
        inputs
    )(calc_data)

    # Plot data
    callback(
        Output(plot_div.plot, 'figure', allow_duplicate=True),
        Input(plot_div.store, 'data'),
        prevent_initial_call='initial_duplicate'
    )(update_plot)

    # Open with WF, PE, and state modals


    # Interaction with WF, PE, and state modals

    # Save data to file
    inputs = [
        Input(options.download_button, 'n_clicks')
    ]
    states = [
        State(plot_div.store, 'data')
    ]
    callback(
        Output(options.download_trigger, 'data'),
        inputs, states, prevent_initial_call=True
    )(download_data)

    # Update plot save format
    outputs = [
        Output(plot_div.plot, 'figure', allow_duplicate=True),
        Output(plot_div.plot, 'config')
    ]
    callback(
        outputs,
        Input(options.image_format_select, 'value'),
        prevent_initial_call=True
    )(update_save_format)

    return


def update_inputs(lin_wn: float, ang_wn: float, fc: float, mu: float,
                  lin_wn_fix: bool, ang_wn_fix: bool, fc_fix: bool,
                  mu_fix: bool) -> tuple[dict[str, float,], list[float], list[bool]]: # noqa
    '''
    Updates input values using checkboxes and current values
    '''

    # Set all text entries as uneditable
    lin_wn_disable = True
    ang_wn_disable = True
    fc_disable = True
    mu_disable = True

    # Make 'fixed' values editable
    if lin_wn_fix:
        lin_wn_disable = False
    if ang_wn_fix:
        ang_wn_disable = False
    if fc_fix:
        fc_disable = False
    if mu_fix:
        mu_disable = False

    if sum([ang_wn_fix, lin_wn_fix, fc_fix, mu_fix]) != 2 or ang_wn_fix and lin_wn_fix: # noqa

        rounded = [
            round(lin_wn, 2), round(ang_wn, 2), round(fc, 2), round(mu, 4)
        ]

        on_off = [
            lin_wn_disable, ang_wn_disable, fc_disable, mu_disable
        ]

        return {'fc': None, 'mu': None}, *rounded, *on_off

    #  Calculate missing parameters
    if ang_wn_fix and not lin_wn_fix:
        omega = ang_wn * LIGHT
        lin_wn = ang_wn / (2 * np.pi)
        if fc_fix and not mu_fix:
            mu = calculate_mu(omega, fc)
        elif not fc_fix and mu_fix:
            fc = calculate_k(omega, mu)
    elif lin_wn_fix and not ang_wn_fix:
        omega = lin_wn * 2 * np.pi * LIGHT
        ang_wn = lin_wn * 2 * np.pi
        if fc_fix and not mu_fix:
            mu = calculate_mu(omega, fc)
        elif not fc_fix and mu_fix:
            fc = calculate_k(omega, mu)
    elif not ang_wn_fix and not lin_wn_fix:
        ang_wn = np.sqrt(fc / (mu * 1.6605E-27)) / LIGHT
        lin_wn = ang_wn / (2 * np.pi)

    rounded = [
        round(lin_wn, 2), round(ang_wn, 2), round(fc, 2), round(mu, 4)
    ]

    on_off = [
        lin_wn_disable, ang_wn_disable, fc_disable, mu_disable
    ]

    return {'fc': fc, 'mu': mu * 1.6605E-27, 'ang_wn': ang_wn, 'lin_wn': lin_wn}, *rounded, *on_off # noqa


def calc_data(vars: dict[str, float], max_n: int):
    '''
    Calculates harmonic potential, states, and wavefunctions

    Parameters
    ----------
    vars: dict[str: float]
        Input data for harmonic oscillator
        keys are 'mu', 'fc', 'ang_wn', 'lin_wn'
    max_n: int
        Maximum number of harmonic states to compute

    Returns
    -------
    data : dict[str: list]
        Keys and values:\n
        'x': list[float] Displacement (x) in metres
        'wf': list[list[float]] Harmonic wavefunction(s) at x
        'states': list[float] Harmonic state energies
        'potential': list[float] Harmonic potential at x
    '''

    if None in vars.values():
        return no_update

    if max_n is None:
        max_n = 5

    # Convert wavenumbers to frequency in units of s^-1
    state_e, harmonic_e, displacement, _ = calc_harmonic_energies(
        vars['fc'],
        vars['mu'],
        max_n
    )

    # Convert to cm-1
    # 1 cm-1 = 1.986 30 x 10-23 J
    state_e /= 1.98630E-23
    harmonic_e /= 1.98630E-23

    wf = [
        harmonic_wf(n, displacement, vars['mu'], vars['ang_wn'] * LIGHT)
        for n in range(max_n + 1)
    ]

    final = {
        'x': (displacement * 10E10).tolist(),
        'wf': wf,
        'states': state_e.tolist(),
        'potential': harmonic_e.tolist()
    }

    return final

# , toggle_pe: bool, toggle_wf: bool,
#                 toggle_states: bool, colour_pbe: str, colour_wf: str,
#                 colour_states: str

def update_plot(data: dict[str, list]) -> Patch:
    '''
    Plots harmonic state energies and wavefunctions, and harmonic potential

    Parameters
    ----------
    data : dict[str: list]
        Keys and values:\n
        'x': list[float] Displacement (x) in metres
        'wf': list[list[float]] Harmonic wavefunction(s) at x
        'states': list[float] Harmonic state energies
        'potential': list[float] Harmonic potential at x
    Returns
    -------
    Patch
        Patch graph figure
    '''
    toggle_pe, toggle_wf, toggle_states = True, True, True

    fig = Patch()

    traces = []

    if toggle_pe:

        # Harmonic potential
        traces.append(
            go.Scatter(
                x=data['x'],
                y=data['potential']
            )
        )

    _states = [
        [da, da]
        for da in data['states']
    ]

    if toggle_states:
        _x = [
            min(data['x']), max(data['x'])
        ]
        # States
        for state in _states:
            traces.append(
                go.Scatter(
                    x=_x,
                    y=state,
                    line={
                        'color': 'black'
                    },
                    mode='lines'
                )
            )

    # Plot harmonic wavefunction
    if toggle_wf:

        # WF
        for n, (wf, state) in enumerate(zip(data['wf'], _states)):
            traces.append(
                go.Scatter(
                    x=data['x'],
                    y=[val * 0.01 + state[0] for val in wf]
                )
            )

    # Find nmax + 1 th state and use this energy as y limit
    upen = 2*data['states'][-1] - data['states'][-2]
    lowen = data['states'][0] - (data['states'][-1] - data['states'][-2]) / 2

    fig['layout']['yaxis']['range'] = [lowen, upen]

    fig['data'] = traces

    return fig


def download_data(_nc: int, data: dict) -> dict:
    '''
    Creates output file for harmonic potential energies, state energies,
    and wavefunctions

    Parameters
    ----------
    _nc: int
        Number of clicks on download button. Used only to trigger callback
    data : dict[str: list]
        Keys and values:\n
        'x': list[float] Displacement (x) in metres
        'wf': list[list[float]] Harmonic wavefunction(s) at x
        'states': list[float] Harmonic state energies
        'potential': list[float] Harmonic potential at x
    Returns
    -------
    dict
        Output dictionary used by dcc.Download
    '''

    if not len(data):
        return no_update

    oc = io.StringIO()

    header = (
        'Vibrational wavefunction data calculated using waveplot.com\n'
        'A tool by Jon Kragskow\n'
    )

    oc.write(header)

    oc.write('\nState energies (cm-1)\n')
    for se in data['states']:
        oc.write('{:.6f}\n'.format(se))

    oc.write(
        '\nDisplacement (A), Harmonic potential (cm-1)\n'
    )
    for di, se in zip(data['x'], data['potential']):
        oc.write('{:.6f}, {:.6f}\n'.format(di * 10E10, se))

    oc.write(
        '\nDisplacement (A), Harmonic Wavefunction for n=0, n=1, ...\n'
    )

    for di, row in zip(data['x'], data['wf']):
        oc.write('{:.8f}, '.format(di * 10E10))
        for state_wf in row:
            oc.write('{:.8f}, '.format(state_wf))
        oc.write('\n')

    output = {
        'content': oc.getvalue(),
        'filename': 'waveplot_harmonic_data.dat'
    }

    return output


def update_save_format(fmt: str):
    '''
    Updates save format of plot.

    Parameters
    ----------
    fmt: str {png, svg, jpeg}
        Image format to use

    Returns
    -------
    Patch
        Patched graph figure
    Patch
        Patched graph config
    '''

    # Figures
    # resetting their layout attr redraws the plot
    # which is necessary because editing the config attr (below) does not...
    fig = Patch()
    fig['layout'] = VIB_LAYOUT

    # Config
    config = Patch()

    # Update plot format in config dict
    config['toImageButtonOptions']['format'] = fmt

    return fig, config
