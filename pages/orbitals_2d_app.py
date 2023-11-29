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

from dash import html, dcc, callback_context, register_page, callback, \
    clientside_callback, ClientsideFunction
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output
import plotly.graph_objs as go
import numpy as np
import io

from .core import orbitals as oc
from .core import utils
from .core import common

ID_PREFIX = 'orb2_'
dash_id = common.dash_id(ID_PREFIX)

PAGE_NAME = '2d Orbitals'
PAGE_PATH = '/orbitals2d'
PAGE_IMAGE = 'assets/Orbitals2.png'
PAGE_DESCRIPTION = 'Interactive atomic orbitals'

register_page(
    __name__,
    order=1,
    path=PAGE_PATH,
    name=PAGE_NAME,
    title=PAGE_NAME,
    image=PAGE_IMAGE,
    description=PAGE_DESCRIPTION
)

prefix = utils.dash_id('orb3')

plot_div = oc.Orb2dPlot(ID_PREFIX)

# Make AC options tab and all callbacks
options_2d = oc.OptionsDiv(ID_PREFIX)
# Connect callbacks for plots and options
oc.assemble_2d_callbacks(
    plot_div, options_2d
)

# Layout of webpage
layout = common.make_layout(plot_div.div, options_2d.div)
